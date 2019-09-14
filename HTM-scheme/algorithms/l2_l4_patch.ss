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
  ;;
  ;; "When in doubt, use brute force" [Ken Thompson]
  #|

A Patch represents an array of cortical columns constructed from layers
("algorithms") with defined interconnections. Inputs are provided as
arguments to a compute procedure; layer outputs can be accessed.
  
L2L4 cortical column:
  L2 is a "column pooler" layer, L4 is a "apical tiebreak temporal memory" layer
  (with ss4L4, ss4L23, and p4 cell populations)
  
                              +---------------<- other L2 cc's active cells
                              v                    (from previous time step)
            +-----------------li--+             
            |          L2        id<--+--------> other L2 cc's lateral inputs  
            +--ff---gc--------ac--+   |            (for next time step)
               ^    ^         v       |              
               |    |         +-------+          ac = active cells
               ^    ^         v                  ai = apical input (to p4 only)
            +--ac---pc--------ai--+              bi = basal input
      +---->bi         L4         |              ff = feedforward input
      |     +-------ff------------+              gc = growth candidates  
      |             ^                            id = internal distal input
      |             |                            li = lateral input (per cc)
      ^             ^                            pc = predicted cells
    location      feature

    |#
  
(library (HTM-scheme HTM-scheme algorithms l2_l4_patch)
                                                                                            ;
(export
  make-patch
  patch-l2s
  patch-l4s
  compute
  reset
  reset-seq
  (rename
    (l4:l4-pop                     l4-pop)
    (l4:get-active-cells           get-active-cells)
    (l4:get-predicted-cells        get-predicted-cells)
    (l4:get-predicted-active-cells get-predicted-active-cells)
    (l4:number-of-basal-segments   number-of-basal-segments)
    (l4:number-of-basal-synapses   number-of-basal-synapses)
    (l4:number-of-apical-segments  number-of-apical-segments)
    (l4:number-of-apical-synapses  number-of-apical-synapses)
    (l2:test:number-of-proximal-segments  get-n-l2-proximal-segs)
    (l2:test:number-of-proximal-synapses  get-n-l2-proximal-syns)
    (l2:test:number-of-distal-segments    get-n-l2-distal-segs)
    (l2:test:number-of-distal-synapses    get-n-l2-distal-syns)
    (l2:test:number-of-lateral-segments   get-n-l2-lateral-segs)
    (l2:test:number-of-lateral-synapses   get-n-l2-lateral-syns)
    ))
                                                                                            ;
(import
  (except (chezscheme) add1 make-list random reset)
          (HTM-scheme HTM-scheme algorithms htm_prelude)
          (HTM-scheme HTM-scheme algorithms htm_concept)
  (prefix (HTM-scheme HTM-scheme algorithms layer4)        l4:)
  (prefix (HTM-scheme HTM-scheme algorithms column_pooler) l2:))
                                                                                            ;                                                                                            ;
(define-record-type patch (fields        ;; Patch
  l2s                                    ;; CColVecOf CP
  l4s)                                   ;; CColVecOf L4
  (protocol
    (lambda (new)
      (lambda (ncc cc ncl4pop nib bis cp-overrides attm-overrides)
        (let* (
            (l2-cell-count (int<- (* 3 cc ncl4pop)))
            (l2s (build-vector ncc (lambda _ 
                    (l2:make-cp (append `([cell-count . ,l2-cell-count]) cp-overrides)))))
            (ais (l2:number-of-cells (vector-ref l2s 0)))
            (l4s (build-vector ncc (lambda _ 
                    (l4:make-l4 cc ncl4pop nib bis ais attm-overrides)))))
          (new l2s l4s))))))
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

(define (lateral-inputs ccx ls prev)     ;; Nat (CCVecOf Layer) (CCVecOf {CellX}) -> { {CellX} }
  ;; produce list of active cell lists for ccs adjacent to ccx from prev timestep
  (let ((layer (vector-ref ls ccx))
        (ncc   (vector-length ls)))
    (vector-fold-left
      (lambda (acc other-layer other-ccx)
        (let ((other-prev-active (vector-ref prev other-ccx)))
          (if (eq? layer other-layer)  acc
            (cond
              [(fx=? ncc 7)                ;; 7 ccs: connect to adjacent cc only
                (if (or (zero? ccx) (zero? other-ccx))
                  (cons other-prev-active acc)
                  (if (fx=? 1 (fxmod (fx- ccx other-ccx) (fx- ncc 1)))
                      (cons other-prev-active acc)
                      acc))]
              [(fx=? ncc 19)               ;; 19 ccs
                (if (memv other-ccx (vector-ref adjacent-ccs ccx))
                    (cons other-prev-active acc)
                    acc)]
              [else
                (cons other-prev-active acc)]))))
      '()
      ls (indexes ls))))
                                                                                            ;
(define (compute patch features locations learn)
  ;; run one timestep of patch
  (let ((l2-prev-actives
          (vector-map (lambda (l2)
              (append                    ;; copy *elements* for next step
                (l2:get-active-cells l2) '()))
            (patch-l2s patch))))
    (threaded-vector-for-each              ;; thread per cortical column
      (lambda (l2 l4 feature location ccx)
        (let-values ([(l4-active-cells l4-predicted-cells bursting-columns)
            (l4:compute l4 
              feature                      ;; feedforward input
              location                     ;; basal input
              (l2:get-active-cells l2)     ;; apical input 
              learn)])
          (when (null? bursting-columns)
            (l2:compute l2
              l4-active-cells              ;; feedforward input
              (lateral-inputs ccx (patch-l2s patch) l2-prev-actives)
              l4-predicted-cells           ;; feedforward growth candidates
              learn
              l4-predicted-cells))))       ;; predicted input
      (patch-l2s patch) (patch-l4s patch) features locations (indexes features))))
                                                                                            ;
(define (reset patch)
  (vector-for-each  l2:reset  (patch-l2s patch))
  (vector-for-each  l4:reset  (patch-l4s patch)))
    
(define (reset-seq patch)
  (vector-for-each  l4:reset-seq  (patch-l4s patch)))
                                                                                            ;
)
