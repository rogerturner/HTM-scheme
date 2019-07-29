#!chezscheme

;; === HTM-scheme L2L4 Patch algorithm Copyright 2019 Roger Turner. ===
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; This program is free software: you can redistribute it and/or modify  ;;
  ;; it under the terms of the GNU Affero Public License version 3 as      ;;
  ;; published by the Free Software Foundation.                            ;;
  ;;                                                                       ;;
  ;; This program is distributed in the hope that it will be useful,       ;;
  ;; but WITHOUT ANY WARRANTY; without even the implied warranty of        ;;
  ;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                  ;;
  ;; See the GNU Affero Public License for more details.                   ;;
  ;;                                                                       ;;
  ;; You should have received a copy of the GNU Affero Public License      ;;
  ;; along with this program.  If not, see http://www.gnu.org/licenses.    ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; This L2L4 patch is part of the untangling_sequences project.
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.
  #|

A Patch represents an array of cortical columns constructed from layers
("algorithms") with defined interconnections. Inputs are provided as
arguments to a compute procedure; layer outputs can be accessed.
  
L2L4 cortical column:

                              +----------+      L2obj is a "column pooler" layer
                              v          |      L4loc and L4seq are "apical tiebreak
            +-----------------li--+      |      temporal memory" layers
            |       L2obj        id<-+   |       
            +--ff---gc--------ac--+  |   +-----<- other L2 cc's active cells            
               ^    ^         v      |              (from previous time step)
         +-----+    |         +------+----------> other L2 cc's lateral inputs     
         |     ^    ^         v      |              (for next time step)
         |  +--ac---pc---+-+--ai--+  |
      +--c->bi   L4loc   |*|      |  |
      |  |  +------------|*|------+  |            ac = active cells
      |  |               |*|         |            ai = apical input
      |  +-----+         |*|  +------+            bi = basal input
      |  |     ^         |*|  v                   ff = feedforward input
      |  |  +--ac--------|*|--ai--+                    (union of L4loc and L4seq ac)
      |  +->bi   L4seq   |*|      |               gc = growth candidates
      |     +------------+-+------+               id = internal distal input
      |                   ^ \                     li = lateral input (per cc)
      |                   |  \                    pc = predicted cells
      ^                   ^   \                   
    location        feature   minicolumn spans L4loc and L4seq

    |#
  
(library (HTM-scheme HTM-scheme algorithms L2L4_patch)
                                                                                            ;
(export
  make-patch
  patch-L2s
  patch-L4s
  compute
  reset
  reset-seq
  (rename
    (l4:L4                         L4)
    (l4:get-active-cells           get-active-cells)
    (l4:get-predicted-cells        get-predicted-cells)
    (l4:get-predicted-active-cells get-predicted-active-cells)
    (l4:get-n-segments-created     get-n-segments-created)
    (l4:get-n-synapses-created     get-n-synapses-created)))
                                                                                            ;
(import
  (except (chezscheme) add1 make-list random reset)
          (HTM-scheme HTM-scheme algorithms htm_prelude)
          (HTM-scheme HTM-scheme algorithms htm_concept)
  (prefix (HTM-scheme HTM-scheme algorithms L4_ss4_a23_attm)  l4:)
  (prefix (HTM-scheme HTM-scheme algorithms column_pooler)    l2:))
                                                                                            ;                                                                                            ;
(define-record-type patch (fields        ;; Patch
  L2s                                    ;; CColVecOf CP
  L4s                                    ;; CColVecOf L4
  (mutable L2-prev-actives)              ;; CColVecOf (listof CellX)
  enable-feedback)                       ;; Boolean
  (protocol
    (lambda (new)
      (lambda (ncc is nib ef cp-overrides loc-overrides seq-overrides)
        (let* ( 
            (L2s (build-vector ncc (lambda _ 
                    (l2:make-cp (append cp-overrides (get-default-cp-params is nib))))))
            (ais (l2:number-of-cells (vector-ref L2s 0)))
            (L4s (build-vector ncc (lambda _ 
                    (l4:make-l4 is nib ais loc-overrides seq-overrides)))))
          (new L2s
               L4s
               (build-vector ncc (lambda _ '()))
               ef))))))
                                                                                            ;
(define adjacent-ccs (vector
  '( 1  2  3  4  5  6)
  '( 7  8  2  0  6 18)                   ;;          17    18     7
  '( 8  9 10  3  0  1)                   ;;
  '( 2 10 11 12  4  0)                   ;;       16     6     1      8
  '( 0  3 12 13 14  5)                   ;;
  '( 6  0  4 14 15 16)                   ;;    15     5     0     2      9
  '(18  1  0  5 16 17)                   ;;
  '( 8  2  1  6 18)                      ;;       14     4     3     10
  '( 9  2  0  1  7)                      ;;
  '(10  3  2  1  8)                      ;;          13    12    11
  '(11  3  0  2  9)
  '(12  4  3  2 10)
  '(13  4  0  3 11)
  '(14  5  4  3 12)
  '(15  5  0  4 13)
  '(16  6  5  4 14)
  '(17  6  0  5 15)
  '(18  1  6  5 16)
  '( 7  1  0  6 17)))

(define (lateral-inputs patch l2 ccx)        ;; -> { {CellX} }
  ;; produce list of active cell lists for adjacent ccs from prev timestep
  (vector-fold-left
    (lambda (lateral-inputs other-l2 other-prev-active other-ccx)
      (let ((nccm1 (fx- (vector-length (patch-L2s patch)) 1)))
        (if (eq? l2 other-l2)  lateral-inputs
          (cond
            [(fx=? nccm1 6)          ;; 7 ccs: connect to adjacent cc only
              (if (or (zero? ccx) (zero? other-ccx))
                (cons other-prev-active lateral-inputs)
                (if (fx=? 1 (fxmod (fx- ccx other-ccx) nccm1))
                    (cons other-prev-active lateral-inputs)
                    lateral-inputs))]
            [(fx=? nccm1 18)         ;; 19 ccs
              (if (memv other-ccx (vector-ref adjacent-ccs ccx))
                  (cons other-prev-active lateral-inputs)
                  lateral-inputs)]
            [else
              (cons other-prev-active lateral-inputs)]))))
    '()
    (patch-L2s patch) (patch-L2-prev-actives patch) (indexes (patch-L2s patch))))
                                                                                            ;
(define compute 
(case-lambda [ (patch features locations learn)
                 (compute patch features locations learn (l4:L4 loc seq)) ]
  [ (patch features locations learn l4pop)
  ;; run one timestep of patch; l4pop enables L4 location/sequence populations
  ;; feature input not predicted by either L4 or TM causes bursting across both;
  ;; if loc & seq then L2 input is union of L4 and TM active cells
  (threaded-vector-for-each                    ;; thread per cortical column
    (lambda (l2 l4 feature location ccx)
      (let ((apical-input  (if (patch-enable-feedback patch)
                                 (l2:get-active-cells l2)  '())))
        (l4:compute l4 feature location apical-input learn l4pop )
        (l2:compute l2
          (l4:get-active-cells l4 l4pop)       ;; feedforward input
          (lateral-inputs patch l2 ccx)
          (l4:get-predicted-cells l4 l4pop)    ;; feedforward growth candidates
          learn
          (l4:get-predicted-cells l4 l4pop)))) ;; predicted input
    (patch-L2s patch) (patch-L4s patch) features locations (indexes features))
  (patch-L2-prev-actives-set! patch            ;; save L2 activity for t+1 lateral input
    (vector-map (lambda (l2)
        (append                                ;; copy *elements* for next step
          (l2:get-active-cells l2) '()))
      (patch-L2s patch))) ] ) )
                                                                                            ;
(define (reset patch)
  (vector-for-each  l2:reset  (patch-L2s patch))
  (vector-for-each  l4:reset  (patch-L4s patch))
  (patch-L2-prev-actives-set! patch
    (build-vector (vector-length (patch-L2s patch)) (lambda _ '()))))
    
(define (reset-seq patch)
  (vector-for-each  l4:reset-seq  (patch-L4s patch)))
                                                                                            ;
(define (get-default-cp-params input-size num-input-bits)
  ;; getDefaultL2Params from l2_l4_inference.py
  (let* ( (sample-size-proximal
            (case num-input-bits
              [(20) 10] [(10) 6] [else (int<- (* num-input-bits .6))]))
          (min-threshold-proximal
            (case num-input-bits
              [(20) 5] [(10) 3] [else (int<- (* sample-size-proximal .6))])))
    `(
      [input-width                   . ,(* input-size 16)]
      [cell-count                    . 4096]
      [sdr-size                      . 40]
      [syn-perm-proximal-inc         . ,(perm 0.1)]
      [syn-perm-proximal-dec         . ,(perm 0.001)]
      [initial-proximal-permanence   . ,(perm 0.6)]
      [min-threshold-proximal        . ,min-threshold-proximal]
      [sample-size-proximal          . ,sample-size-proximal]
      [connected-permanence-proximal . ,(perm 0.5)]
      [syn-perm-distal-inc           . ,(perm 0.1)]
      [syn-perm-distal-dec           . ,(perm 0.001)]
      [initial-distal-permanence     . ,(perm 0.41)]
      [activation-threshold-distal   . 13]
      [sample-size-distal            . 20]
      [connected-permanence-distal   . ,(perm 0.5)])))

)