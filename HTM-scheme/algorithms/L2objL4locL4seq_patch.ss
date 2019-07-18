#!chezscheme

;; === HTM-scheme L2objL4locL4seq Patch algorithm Copyright 2019 Roger Turner. ===
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

  ;; This L2objL4locL4seq patch is part of the untangling_sequences project.
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.
  #|

A Patch represents an array of cortical columns constructed from layers
("algorithms") with defined interconnections. Inputs are provided as
arguments to a compute procedure; layer outputs can be accessed.
  
L2objL4locL4seq cortical column:

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
  
(library (HTM-scheme HTM-scheme algorithms L2objL4locL4seq_patch)
                                                                                            ;
(export
  make-patch
  patch-L2s
  patch-L4s
  patch-TMs
  compute
  reset
  L4)
                                                                                            ;
(import
  (except (chezscheme) add1 make-list random reset)
          (HTM-scheme HTM-scheme algorithms htm_prelude)
          (HTM-scheme HTM-scheme algorithms htm_concept)
  (prefix (HTM-scheme HTM-scheme algorithms apical_tiebreak_temporal_memory) attm:)
  (prefix (HTM-scheme HTM-scheme algorithms column_pooler)                   cp:))
                                                                                            ;
(define-enumeration L4-populations (loc seq) L4)
                                                                                            ;
(define-record-type patch (fields        ;; Patch
  L2s                                    ;; CColVecOf CP
  L4s                                    ;; CColVecOf ATTM
  TMs                                    ;; CColVecOf ATTM
  (mutable L2-prev-actives)              ;; CColVecOf (listof CellX)
  enable-feedback)                       ;; Boolean
  (protocol
    (lambda (new)
      (lambda (ncc is nib ef L2-overrides L4-overrides TM-overrides)
        (let* ( (L2s (build-vector ncc (lambda _ 
                        (cp:make-cp (append L2-overrides (get-default-L2-params is nib))))))
                (L4s (build-vector ncc (lambda _ 
                        (attm:make-tm (append
                          `([basal-input-size  . ,is]
                            [apical-input-size . ,(cp:number-of-cells (vector-ref L2s 0))])
                          L4-overrides (get-default-L4-params is nib))))))
                (TMs (build-vector ncc (lambda _ 
                        (attm:make-tm (append
                          `([basal-input-size  . ,is]
                            [apical-input-size . ,(cp:number-of-cells (vector-ref L2s 0))])
                          TM-overrides (get-default-TM-params is nib)))))))
          (new L2s L4s TMs (build-vector ncc (lambda _ '())) ef))))))
                                                                                            ;
(define adjacent-ccs (vector
  '( 1  2  3  4  5  6)                   ;;                18    
  '( 7  8  2  0  6 18)                   ;;          17           7
  '( 8  9 10  3  0  1)                   ;;
  '( 2 10 11 12  4  0)                   ;;      16      6     1      8
  '( 0  3 12 13 14  5)                   ;;
  '( 6  0  4 14 15 16)                   ;;     15    5     0     2     9
  '(18  1  0  5 16 17)                   ;;
  '( 8  2  1  6 18)                      ;;       14     4     3     10
  '( 9  2  0  1  7)                      ;;
  '(10  3  2  1  8)                      ;;          13          11
  '(11  3  0  2  9)                      ;;                12
  '(12  4  3  2 10)
  '(13  4  0  3 11)
  '(14  5  4  3 12)
  '(15  5  0  4 13)
  '(16  6  5  4 14)
  '(17  6  0  5 15)
  '(18  1  6  5 16)
  '( 7  1  0  6 17)))

(define compute 
(case-lambda [ (patch features locations learn)
                 (compute patch features locations learn (L4 loc seq)) ]
  [ (patch features locations learn L4pop)
  ;; run one timestep of patch; L4pop enables L4 location/sequence populations
  ;; feature input not predicted by either L4 or TM causes bursting across both;
  ;; if loc & seq then L2 input is union of L4 and TM active cells
  (threaded-vector-for-each                    ;; thread per cortical column
    (lambda (L4 L2 TM feature location ccx)
      (let* ( (loc                     (enum-set-member? 'loc L4pop))
              (seq                     (enum-set-member? 'seq L4pop))
              (basal-input             (if seq (attm:get-active-cells TM)  '()))
              (basal-growth-candidates (if seq (attm:get-winner-cells TM)  '()))
              (apical-input            (if (patch-enable-feedback patch)
                                           (cp:get-active-cells L2)  '())))
        (assert (or loc seq))
        (when loc (attm:depolarize-cells L4 location    apical-input learn))
        (when seq (attm:depolarize-cells TM basal-input apical-input learn))
        (let* (
            (L4-predicted-cells  (if loc (attm:get-predicted-cells L4)  '()))
            (L4-predicted-cols   (attm:cols-from-cells L4 L4-predicted-cells))
            (TM-predicted-cells  (if loc (attm:get-predicted-cells TM)  '()))
            (TM-predicted-cols   (attm:cols-from-cells TM TM-predicted-cells))
            (bursting-columns   
              (setdiff1d feature (cond
                  [(and loc seq)
                    (unique! fx=? (union1d L4-predicted-cols TM-predicted-cols))]
                  [loc  L4-predicted-cols]
                  [seq  TM-predicted-cols])))
            (predicted-cells  (cond
                  [(and loc seq)
                    (unique! fx=? (union1d L4-predicted-cells TM-predicted-cells))]
                  [loc  L4-predicted-cells]
                  [seq  TM-predicted-cells])))
#;>         (define (lateral-input)        ;; -> { {CellX} }
              ;; produce list of active cell lists for adjacent ccs from prev timestep
              (vector-fold-left
                (lambda (lateral-inputs other-L2 other-prev-active other-ccx)
                  (let ((nccm1 (fx- (vector-length (patch-L2s patch)) 1)))
                    (if (eq? L2 other-L2)  lateral-inputs
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
                (patch-L2s patch) (patch-L2-prev-actives patch) (indexes features)))
                                                                                            ;
          (when loc (attm:activate-cells L4 feature location apical-input
                      location apical-input learn bursting-columns))
          (when seq (attm:activate-cells TM feature basal-input apical-input
                      basal-growth-candidates apical-input learn bursting-columns))
          (cp:compute L2
            (cond                              ;; feedforward input
              [(and loc seq)
                (unique! fx=? (union1d (attm:get-active-cells L4) (attm:get-active-cells TM)))]
              [loc (attm:get-active-cells L4)]
              [seq (attm:get-active-cells TM)])
            (lateral-input)
            predicted-cells                    ;; feedforward growth candidates
            learn
            predicted-cells))))                ;; predicted input
    (patch-L4s patch) (patch-L2s patch) (patch-TMs patch) features locations (indexes features))
  (patch-L2-prev-actives-set! patch            ;; save L2 activity for t+1 lateral input
    (vector-map (lambda (L2)
        (append                                ;; copy *elements* for next step
          (cp:get-active-cells L2) '()))
      (patch-L2s patch))) ] ) )
                                                                                            ;
(define (reset patch)
  (vector-for-each
    (lambda (L2 L4 TM)
      (cp:reset   L2)
      (attm:reset L4)
      (attm:reset TM))
    (patch-L2s patch) (patch-L4s patch) (patch-TMs patch))
  (patch-L2-prev-actives-set! patch
    (build-vector (vector-length (patch-L2s patch)) (lambda _ '()))))
                                                                                            ;
(define (get-default-L2-params input-size num-input-bits)
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
                                                                                            ;
(define (get-default-L4-params input-size num-input-bits)
  ;; getDefaultL4Params from l2_l4_inference.py
  (let* ( (sample-size (int<- (* 1.5 num-input-bits)))
          (activation-threshold 
            (case num-input-bits
              [(20) 13] [(10) 8] [else (int<- (* num-input-bits .6))]))
          (min-threshold
            (case num-input-bits
              [(20) 13] [(10) 8] [else activation-threshold])))
    `(
      [column-count                       . ,input-size]
      [cells-per-column                   . 16]
      [initial-permanence                 . ,(perm 0.51)]
      [connected-permanence               . ,(perm 0.6)]
      [permanence-increment               . ,(perm 0.1)]
      [permanence-decrement               . ,(perm 0.02)]
      [min-threshold                      . ,min-threshold]
      [basal-predicted-segment-decrement  . ,(perm 0.0)]
      [apical-predicted-segment-decrement . ,(perm 0.0)]
      [activation-threshold               . ,activation-threshold]
      [reduced-basal-threshold            . ,(int<- (* activation-threshold 0.6))]
      [sample-size                        . ,sample-size])))
                                                                                            ;
(define (get-default-TM-params input-size num-input-bits)
  ;; getDefaultTMParams from combined_sequence_experiment.py
  (let* ( (sample-size (int<- (* 1.5 num-input-bits)))
          (activation-threshold 
            (case num-input-bits
              [(20) 18] [(10) 8] [else (int<- (* num-input-bits .6))]))
          (min-threshold
            (case num-input-bits
              [(20) 18] [(10) 8] [else activation-threshold])))
    `(
      [column-count                       . ,input-size]
      [cells-per-column                   . 16]
      [initial-permanence                 . ,(perm 0.41)]
      [connected-permanence               . ,(perm 0.6)]
      [permanence-increment               . ,(perm 0.1)]
      [permanence-decrement               . ,(perm 0.03)]
      [min-threshold                      . ,min-threshold]
      [basal-predicted-segment-decrement  . ,(perm 0.003)]
      [apical-predicted-segment-decrement . ,(perm 0.0)]
      [reduced-basal-threshold            . ,(int<- (* activation-threshold 0.6))]
      [activation-threshold               . ,activation-threshold]
      [sample-size                        . ,sample-size])))
                                                                                            ;
)