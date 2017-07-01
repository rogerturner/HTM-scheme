#!r6rs

;; ========= HTM-scheme Temporal Memory Copyright 2017 Roger Turner. =========
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Based on temporal_memory.py which is part of the Numenta Platform for ;;
  ;; Intelligent Computing (NuPIC) Copyright (C) 2014-2016, Numenta, Inc.  ;;
  ;;                                                                       ;;
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

  ;; Translated from NuPIC temporal_memory.py, see comments there for more info.
  ;; Numbered comments in the "Temporal Memory Algorithm" section echo, where
  ;; possible, pseudocode in Numenta BaMI Temporal-Memory-Algorithm-Details.pdf 
  ;; revision 0.5
  ;; Use a "Fold All" view (in eg Atom) for a source overview.

(library (libraries lib-tm)

  (export
    temporal-memory
    tm-perm
    tm-compute
    tm-reset
    tm-get-active-cols
    tm-get-predictive-cols)
    
  (import 
    (rnrs)                    ;; use (except (chezscheme) add1 make-list random) for Chez load-program
    (libraries htm-prelude))
    
  #| comment out for Chez Scheme (optional) |#
  (define fxvector        vector)
  (define make-fxvector   make-vector)
  (define fxvector-ref    vector-ref)
  (define fxvector-set!   vector-set!)
  (define fxvector-length vector-length)
  (define list->fxvector  list->vector)
  (define fxvector->list  vector->list)
  #| |#
          
  (define synapses:       fxvector)
  (define make-synapses   make-fxvector)
  (define synapses-ref    fxvector-ref)
  (define synapses-set!   fxvector-set!)
  (define synapses-length fxvector-length)
  (define list->synapses  list->fxvector)
  (define synapses->list  fxvector->list)

;; === Temporal Memory Types ===
                                                                                            ;
;; X, Y, Z      = type parameters (arbitrary types in function type specification etc)
;; X Y -> Z     = function with argument types X, Y and result type Z
;; {X}          = abbreviation for (listof X)
;; Nat          = natural number (including zero) (Scheme Fixnum or exact Integer)
;; Fix10k       = Fixnum interpreted as number with 4 decimal places, ie 10000 = 1.0
;; Permanence   = Fix10k [0 - MAXPERM], interpreted as 0.0 - 0.9999
;; Synapse      = Fixnum, interpreted as CellX * 10000 + Permanence 
;; Synapses     = (vectorof Synapse) [fxvector in Chez Scheme]
;; Segment      = Record: CellX, FlatX, and Synapses
;; FlatX        = Nat, segment sequence number
;; Cell         = (listof Segment)
;; CellX        = Nat [0 - (div (expt 2 (fixnum-width)) 10000)] cell index
;; CellVecOf    = Vector indexed by CellX
;; ColX         = Nat [0 - MAXCOL], column index of cell
;; TM           = Record: tm parameters, cells

;; === Parameters, Data, Convenience Functions ===
                                                                                            ;
(define-record-type tm                   ;; TM
  (fields
    column-dimensions                    ;; (listof Nat)
    cells-per-column                     ;; Nat
    activation-threshold                 ;; Nat
    initial-permanence                   ;; Permanence
    connected-permanence                 ;; Permanence
    min-threshold                        ;; Nat
    max-new-synapse-count                ;; Nat
    permanence-increment                 ;; Permanence
    permanence-decrement                 ;; Permanence
    predicted-segment-decrement          ;; Permanence
    max-segments-per-cell                ;; Nat
    max-synapses-per-segment             ;; Nat
    (mutable next-flatx)                 ;; Nat: next available index in seg-register
    (mutable free-flatx)                 ;; (listof Nat): indices available for re-use
    (mutable active-cells)               ;; (listof CellX)
    (mutable winner-cells)               ;; (listof CellX)
    (mutable active-segments)            ;; (listof Segment)
    (mutable matching-segments)          ;; (listof Segment)
    (mutable iteration)                  ;; Nat
    num-columns                          ;; Nat
    num-cells                            ;; Nat
    cells                                ;; (CellVecOf {Segment})
    pre-index                            ;; (CellVecOf {FlatX})
    (mutable seg-register)               ;; (vectorof Segment)
    (mutable num-active-pot-syns-for-seg);; (hash FlatX -> Nat)
    ))
                                                                                            ;
(define-record-type seg                  ;; Segment
  (fields
    cellx                                ;; Nat: index of cell that this is a segment of
    flatx                                ;; Nat: index of this in the register of segments
    (mutable last-used)                  ;; Nat: iteration number last time segment active
    (mutable synapses)                   ;; (vectorof Synapse)
    ))
                                                                                            ;
(define %max-perm% 9999)                 ;; Fixnum
(define %min-perm% 0)                    ;; Fixnum
                                                                                            ;
(define (clip-p perm)                    ;; Permanence -> Permanence
  (fxmax %min-perm% (fxmin %max-perm% perm)))
                                                                                            ;
(define (tm-perm x)                      ;; Number[0.0-.9999] -> Permanence
  (clip-p (int<- (* x x10k))))
                                                                                            ;
(define tm-defaults                      ;; (listof KWarg)
  `(
    [activation-threshold        . 13]
    [initial-permanence          . ,(tm-perm 0.21)]
    [connected-permanence        . ,(tm-perm 0.5)]
    [min-threshold               . 10]
    [max-new-synapse-count       . 20]
    [permanence-increment        . ,(tm-perm 0.1)]
    [permanence-decrement        . ,(tm-perm 0.1)]
    [predicted-segment-decrement . ,(tm-perm 0.0004)]
    [max-segments-per-cell       . 255]
    [max-synapses-per-segment    . 255]
    [next-flatx                  . 0]
    [free-flatx                  . ()]
    [active-cells                . ()]
    [winner-cells                . ()]
    [active-segments             . ()]
    [matching-segments           . ()]
    [iteration                   . 0]
    ))
                                                                                            ;
(define (temporal-memory cd cpc . args)  ;; (listof Nat) Nat (listof KWarg) -> TM
  ;; produce temporal memory instance with defaults for parameters not specified
  (let* ((num-cells (* cpc (apply * cd)))
         (tm (apply make-tm (append (list cd cpc) 
                           (key-word-args args tm-defaults)
                           (list (apply * cd)                 ;; num-columns
                                 num-cells                    ;; num-cells
                                 (make-vector num-cells '())  ;; cells
                                 (make-vector num-cells '())  ;; pre-index
                                 #f                           ;; seg-register
                                 #f                           ;; num-active-pot-syns-for-seg
                                 )))))
    (tm-seg-register-set! tm (make-vector (* num-cells (tm-max-segments-per-cell tm)) #f))
    tm))
    
  [define TM22 [temporal-memory '(2) 2 '[max-segments-per-cell . 2]] ]
  [expect ([vector-length [tm-seg-register TM22]] 8 )]
                                                                                            ;
;; --- Segments ---
                                                                                            ;
(define (create-segment tm cellx)        ;; TM CellX -> Seg
  ;; produce a new segment, updating the register of segments
  (let ((segments (vector-ref (tm-cells tm) cellx)))
    (when (fx>=? (length segments) (tm-max-segments-per-cell tm))
      (let loop ((segs segments) (oldest-seg #f) (last-used (greatest-fixnum)))
        (cond [ (null? segs)
                (vector-set! (tm-cells tm) cellx (remq oldest-seg segments))
                (destroy-segment tm (car oldest-seg)) ]
              [ (fx<? (seg-last-used (car segs)) last-used)
                (loop (cdr segs) segs (seg-last-used (car segs))) ]
              [ else (loop (cdr segs) oldest-seg last-used) ])))
    (let* ( (new-flatx (if (null? (tm-free-flatx tm))
                            (tm-next-flatx tm)
                            (car (tm-free-flatx tm))))
            (segment (make-seg cellx new-flatx (tm-iteration tm) (make-synapses 0))))
      (vector-set! (tm-seg-register tm) new-flatx segment)
      (if (null? (tm-free-flatx tm))
        (tm-next-flatx-set! tm (add1 new-flatx))
        (tm-free-flatx-set! tm (cdr (tm-free-flatx tm))))
      (vector-set! (tm-cells tm) cellx (cons segment segments))
      segment)))
    
  [define SEG00 [create-segment TM22 0] ]
  [expect ([tm-next-flatx TM22] 1 )]
                                                                                            ;
(define (destroy-segment tm segment)     ;; TM Segment ->
  ;; make segment's flatx available for re-use
  (let* ( (flatx (seg-flatx segment))
          (null-segment (make-seg -1 flatx 0 (make-synapses 0))))
    (vector-set! (tm-seg-register tm) flatx null-segment)
    (tm-free-flatx-set! tm (cons flatx (tm-free-flatx tm)))))
    
  [expect ([destroy-segment TM22 SEG00] #f )
          ([tm-free-flatx TM22] '(0) )
          ([seg-cellx (vector-ref [tm-seg-register TM22] 0)] -1 )
          ([seg-flatx [create-segment TM22 1]] 0 )
          ([tm-free-flatx TM22] '() )]
                                                                                            ;
(define (cellx->colx tm cellx)           ;; TM CellX -> ColX
  (fxdiv cellx (tm-cells-per-column tm)))
                                                                                            ;
(define (seg-colx tm segment)            ;; TM Segment -> ColX
  (cellx->colx tm (seg-cellx segment)))
                                                                                            ;
(define (segs->colx tm segments)         ;; TM (listof Segment) -> ColX
  ;; produce column index for first segment in list, or terminating value if list null
  (if (null? segments) (tm-num-columns tm)
      (cellx->colx tm (seg-cellx (car segments)))))
                                                                                            ;
(define (add-to-pre-index tm prex seg)   ;; TM CellX Segment ->
  ;; add to the list of segment flatxs with synapses from the pre-synaptic cell
  (vector-set! (tm-pre-index tm) prex 
    (cons (seg-flatx seg) (vector-ref (tm-pre-index tm) prex))))
                                                                                            ;
;; --- Permanences and Synapses ---
                                                                                            ;
(define (make-synapse prex perm)         ;; PreX Permanence -> Synapse
  (+ (* prex x10k) (least-fixnum) perm))
                                                                                            ;
(define (prex synapse)                   ;; Synapse -> PreX
  (div (- synapse (least-fixnum)) x10k))
                                                                                            ;
(define (perm synapse)                   ;; Synapse -> Permanence
  (mod (- synapse (least-fixnum)) x10k))
                                                                                            ;
(define (build-synapses n proc)          ;; Nat (Nat -> Synapse) -> Synapses
  (let ((synapses (make-synapses n)))
    (do ((i 0 (add1 i))) ((fx=? i n) synapses)
      (synapses-set! synapses i (proc i)))))
                                                                                            ;
(define (memv-prex? cellx synapses)      ;; CellX (vectorof Synapse) -> Bool
  ;; produce whether cellx equals any pre-synaptic-cell index of synapses
  (let loop ((sx (fx- (synapses-length synapses) 1)))
    (cond [ (negative? sx) #f ]
          [ (fx=? cellx (prex (synapses-ref synapses sx))) #t ]
          [ else (loop (fx- sx 1)) ])))
          
  [expect ([memv-prex? 1 (synapses: [make-synapse 0 2] [make-synapse 1 2]) ] #t )
          ([memv-prex? 0 (synapses: [make-synapse 1 0] [make-synapse 1 2]) ] #f )]
                                                                                            ;
(define (colxs->colx tm colxs)           ;; TM (listof ColX) -> ColX
  ;; produce first column index from list, or terminating value if list null
  (if (null? colxs) (tm-num-columns tm) (car colxs)))
                                                                                            ;
(define (prune-synapses synapses omits)  ;; Synapses (listof Nat) -> Synapses
  ;; omit elements indexed by omits (which is sorted) from synapses
  (let* ( (reslen (fx- (synapses-length synapses) (length omits)))
          (result (make-synapses reslen)))
    (let loop ((rx 0) (sx 0) (omits omits))
      (cond [ (fx=? rx reslen) result ]
            [ (if (null? omits) #f
                (fx=? sx (car omits))) (loop rx (add1 sx) (cdr omits)) ]
            [ else
                (synapses-set! result rx (synapses-ref synapses sx))
                (loop (add1 rx) (add1 sx) omits) ]))))
              
  #;(let* ((s make-synapse) (s0 (s 0 0)) (s1 (s 1 1)) (s2 (s 2 2)) (s3 (s 3 3)))
    [expect ([prune-synapses (synapses: s0 s1 s2 s3) '(0 2)    ] (synapses: s1 s3)       )
            ([prune-synapses (synapses: s0 s1 s2 s3) '()       ] (synapses: s0 s1 s2 s3) )
            ([prune-synapses (synapses: s0 s1 s2 s3) '(0 1 2 3)] (synapses: )            )])
                                                                                            ;
(define (destroy-min-permanence-synapses ;; TM Seg Nat (listof CellX) -> 
          tm segment n-destroy exclude-cells)
  ;; drop n-destroy lowest-permanence synapses, but retaining any on exclude list
  (let-values ([(keepers destroy-candidates)
                (partition (lambda (x) (memv (prex x) exclude-cells))
                           (synapses->list (seg-synapses segment)))])
    (let* ( (destroy-candidates 
              (list-sort (lambda (x y) (fx<? (perm x) (perm y))) destroy-candidates))
            (survivors (list-tail destroy-candidates (fxmin n-destroy (length destroy-candidates)))))
      (seg-synapses-set! segment (list->synapses (append survivors keepers))))))
  
  #;(let ((s make-synapse))
    (vector-set! (tm-cells TM22) 0 (make-seg 0 0 0 (synapses: (s 1 2) (s 3 1) (s 4 1) )))    
    [destroy-min-permanence-synapses TM22 (vector-ref (tm-cells TM22) 0) 1 '(3)]
    [expect ([seg-synapses (vector-ref (tm-cells TM22) 0)] (synapses: (s 1 2) (s 3 1)) )])

;; === Temporal Memory Algorithm ===
                                                                                            ;
;; --- Evaluate the active columns against predictions. Choose a set of active cells. ---
                                                                                            ;
(define (activate-cells tm actcols learn);; TM (listof ColumnX) Boolean ->
  ;; step through columns in merge of actcols/actsegs/matsegs lists, building active/winner lists; 
  ;; the activate/burst/punish subfunctions process one column and return list for next column
  ;; actcells/wincells are current; TM active/winner-cells are previous lists until end
  (let loop ( (actcols actcols)                                                ;; 1. for column in columns
              (actsegs (tm-active-segments tm)) 
              (matsegs (tm-matching-segments tm)) 
              (actcells '()) 
              (wincells '()))
    (let* ( (actcols-colx (colxs->colx tm actcols))
            (actsegs-colx (segs->colx  tm actsegs))
            (matsegs-colx (segs->colx  tm matsegs))
            (column (fxmin actcols-colx actsegs-colx matsegs-colx)))
      (if (fx<? column (tm-num-columns tm))          ;; more cols?
        (if (fx=? column actcols-colx)                                         ;; 2. if column in activeColumns(t) then
            (if (fx=? column actsegs-colx)                                     ;; 3.   if count(segmentsForColumn(column, activeSegments(t-1))) > 0 then
                (let-values ([(nextactseg actcells wincells)                   ;; 4.     activatePredictedColumn(column)
                    (activate-predicted-column tm column actsegs actcells wincells learn)])
                  (loop actcols nextactseg matsegs actcells wincells))
                (let-values ([(nextmatseg actcells wincells)                   ;; 6.   else burstColumn(column)
                    (burst-column tm column matsegs actcells wincells learn)])
                  (loop (cdr actcols) actsegs nextmatseg actcells wincells)))
            (if (fx=? column matsegs-colx)                                     ;; 8.  else if count(segmentsForColumn(column, matchingSegments(t-1))) > 0 then
                (let ((nextmatseg (if learn                                    ;; 50.   if LEARNING_ENABLED
                                    (punish-predicted-column tm column matsegs);; 9.      punishPredictedColumn(column)
                                    (skip-col tm column matsegs))))
                  (loop (if (null? actcols) actcols (cdr actcols))
                        (skip-col tm column actsegs)
                        nextmatseg actcells wincells))
                (loop (if (null? actcols) actcols (cdr actcols))
                      (skip-col tm column actsegs)
                      matsegs actcells wincells)))
        ;; all active cols and segs handled: set prev active/winner cells for next iteration
        (begin
          (tm-active-cells-set! tm actcells)
          (tm-winner-cells-set! tm wincells))))))
                                                                                            ;
(define (activate-predicted-column tm    ;; TM ColX {Seg} {CellX} {CellX} Bool -> {Seg} {CellX} {CellX}
          column active-segs actcells wincells learn)                          ;; 10. function activatePredictedColumn(column)
  ;; use active-segs for column to add to active/winner cells lists, 
  ;; step to active-seg for next column
  (let loop ( (segments active-segs) 
              (actcells actcells) 
              (wincells wincells) 
              (prev-cell (tm-num-cells tm)))
    (if (fx=? (segs->colx tm segments) column)                                 ;; 11. for segment in segmentsForColumn(column, activeSegments(t-1))
        (let ((segment (car segments)))
          (when learn                                                          ;; 15. if LEARNING_ENABLED
            (adapt-segment tm (car segments) (tm-active-cells tm)
                (tm-permanence-increment tm) (tm-permanence-decrement tm))
            (let ((n-grow-desired                                              ;; 22. newSynapseCount = (SYNAPSE_SAMPLE_SIZE -
                    (fx- (tm-max-new-synapse-count tm)
                         (hashtable-ref (tm-num-active-pot-syns-for-seg tm)    ;; 23.                    numActivePotentialSynapses(t-1, segment))
                                        (seg-flatx segment) 0))))
              (when (fxpositive? n-grow-desired)
                (grow-synapses tm segment n-grow-desired                       ;; 24. growSynapses(segment, newSynapseCount)
                               (tm-winner-cells tm)))))
          (let ((this-cell (seg-cellx segment)))
            (if (fx=? this-cell prev-cell)
                (loop (cdr segments) actcells wincells prev-cell)
                (loop (cdr segments) 
                      (cons this-cell actcells)                                ;; 12. activeCells(t).add(segment.cell)
                      (cons this-cell wincells)                                ;; 13. winnerCells(t).add(segment.cell)
                      this-cell))))
        ;; end of segments for this column, return next segment and updated cell lists
        (values segments actcells wincells))))
                                                                                            ;
(define (adapt-segment tm                ;; TM {Seg} (listof CellX) Permanence Permanence ->
          segment prev-active-cells permanence-increment permanence-decrement)
  ;; update synapses on segment: + if pre cell active, - otherwise, destroy if 0
  (let ((synapses (seg-synapses segment)))
    (let loop ((sx (fx- (synapses-length synapses) 1)) (synapses-to-destroy '())) ;; 16. for synapse in segment.synapses
      (if (negative? sx)
          (cond [ (null? synapses-to-destroy) ]
                [ (fx=? (length synapses-to-destroy) (synapses-length synapses))
                  (vector-set! (tm-cells tm) (seg-cellx segment)
                    (remq segment (vector-ref (tm-cells tm) (seg-cellx segment))))
                  (destroy-segment tm segment) ]
                [ else (seg-synapses-set! segment 
                        (prune-synapses (seg-synapses segment) synapses-to-destroy)) ])
          (let* ( (synapse (synapses-ref synapses sx))
                  (permanence
                    (if (memv (prex synapse) prev-active-cells)                ;; 17. if synapse.presynapticCell in activeCells(t-1) then
                      (clip-p (fx+ (perm synapse) permanence-increment))       ;; 18.   synapse.permanence += PERMANENCE_INCREMENT
                      (clip-p (fx- (perm synapse) permanence-decrement)))))    ;; 20. else synapse.permanence -= PERMANENCE_DECREMENT
            (if (zero? permanence)
                (loop (fx- sx 1) (cons sx synapses-to-destroy))  ;; build s-t-d sorted
                (begin
                  (synapses-set! synapses sx                                   ;; 18.
                                 (make-synapse (prex synapse) permanence))     ;; 20.
                  (loop (fx- sx 1) synapses-to-destroy))))))))
                                                                                            ;
(define (burst-column tm                 ;; TM ColX {Seg} {CellX} {CellX} Bool -> {Seg} {CellX} {CellX}
          column column-matching-segments actcells wincells learn)             ;; 25. function burstColumn(column)
  ;; add to active/winner cell lists, step to next matching-seg
  (let* ( (start (fx* column (tm-cells-per-column tm)))
          (cells-for-column (build-list (tm-cells-per-column tm)               ;; 26. for cell in column.cells
                                        (lambda (x) (fx+ x start)))))
    (if (fx=? column (segs->colx tm column-matching-segments))                 ;; 29. if segmentsForColumn(column, matchingSegments(t-1)).length > 0
        (let-values ([(nextmatseg best-matching-segment)                       ;; 30.   learningSegment = bestMatchingSegment(column)
                      (best-matching-segment tm column column-matching-segments)])
          (when learn                                                          ;; 39. if LEARNING_ENABLED
            (adapt-segment tm (car best-matching-segment) (tm-active-cells tm) ;; 40-44. [use adapt-segment to adjust permanences]
                (tm-permanence-increment tm) (tm-permanence-decrement tm))
            (let ((n-grow-desired
                    (fx- (tm-max-new-synapse-count tm)                         ;; 46. newSynapseCount = (SAMPLE_SIZE -
                         (hashtable-ref (tm-num-active-pot-syns-for-seg tm)    ;; 47.   numActivePotentialSynapses(t-1, learningSegment))
                                        (seg-flatx (car best-matching-segment)) 0))))
              (when (fxpositive? n-grow-desired)                               ;; 48. growSynapses(learningSegment, newSynapseCount)
                (grow-synapses tm (car best-matching-segment) n-grow-desired (tm-winner-cells tm)))))
          (values nextmatseg 
                  (append cells-for-column actcells)                           ;; 27. activeCells(t).add(cell)
                  (cons (seg-cellx (car best-matching-segment)) wincells)))    ;; 37/31. winnerCells(t).add(learningSegment.cell)

        (let ((winner-cell (least-used-cell tm cells-for-column)))             ;; 33. else winnerCell = leastUsedCell(column)
          (when learn                                                          ;; 34.   if LEARNING_ENABLED
            (let ((n-grow-exact (fxmin (tm-max-new-synapse-count tm)
                                       (length (tm-winner-cells tm)))))
              (when (fxpositive? n-grow-exact)
                (let ((lseg (create-segment tm winner-cell)))                  ;; 35.   learningSegment = growNewSegment(winnerCell)
                  (grow-synapses tm lseg n-grow-exact (tm-winner-cells tm))))));; 48. growSynapses(learningSegment, newSynapseCount)
          (values column-matching-segments 
                  (append cells-for-column actcells)                           ;; 27. activeCells(t).add(cell)
                  (cons winner-cell wincells))))))                             ;; 37. winnerCells(t).add(winnerCell)
                                                                                            ;
(define (punish-predicted-column tm      ;; TM ColX (listof Seg) -> (listof Seg)
          column column-matching-segments)                                     ;; 49. function punishPredictedColumn(column)
  ;; weaken synapses in segments of column; step to next column's segments
  (let loop ((segs column-matching-segments))                                  ;; 51. for segment in segmentsForColumn(column, matchingSegments(t-1))
    (if (fx=? column (segs->colx tm segs))
      (begin
        (when (fxpositive? (tm-predicted-segment-decrement tm))                ;; 52-54. [use adapt-segment to decrement permanences]
          (adapt-segment tm (car segs) (tm-active-cells tm) 
                         (- (tm-predicted-segment-decrement tm)) 0))
        (loop (cdr segs)))
      segs)))
                                                                                            ;
;; --- Activate a set of dendrite segments. ---
                                                                                            ;
(define (compute-activity tm             ;; TM (listof CellX) -> (hash Flatx -> Nat) (hash Flatx -> Nat)
          active-presynaptic-cells)
  ;; produce hashtables mapping indices of segments with pre synapses to count of active cells
  (let ((napsfs (make-eqv-hashtable))  ;; "num-active-potential-synapses-for-segment"
        (nacsfs (make-eqv-hashtable))  ;; "num-active-connected-synapses-for-segment"
        (threshold (tm-connected-permanence tm)))
    (for-each                                                                  ;; 59. if synapse.presynapticCell in activeCells(t) then
      (lambda (cellx)
        (for-each                                                              ;; 55. for segment in segments
          (lambda (flatx)
            (let ((synapses (seg-synapses (vector-ref (tm-seg-register tm) flatx))))
              (do ((sx 0 (add1 sx))) ((fx=? sx (synapses-length synapses)))    ;; 58. for synapse in segment.synapses
                (let ((synapse (synapses-ref synapses sx)))
                  (when (fx=? cellx (prex synapse))                            ;; 63. if synapse.permanence ≥ 0 then
                    (hashtable-update! napsfs flatx add1 0)                    ;; 64.   numActivePotential += 1
                      (when (fx>=? (perm synapse) threshold)                   ;; 60. if synapse.permanence ≥ CONNECTED_PERMANENCE then
                        (hashtable-update! nacsfs flatx add1 0)))))))          ;; 61.   numActiveConnected += 1
          (vector-ref (tm-pre-index tm) cellx)))
      active-presynaptic-cells)
    (values nacsfs napsfs)))
                                                                                            ;
(define (activate-dendrites tm learn)    ;; TM Bool ->
  ;; set active/matching segments lists (sorted by column); set last-used iteration in active segments
  ;; use compute-activity to locate and count segments with synapses for active cells
  (let*-values (
    [(num-active-connected num-active-potential) (compute-activity tm (tm-active-cells tm))]
    [(active-flatxs   active-counts)   (hashtable-entries num-active-connected)]
    [(matching-flatxs matching-counts) (hashtable-entries num-active-potential)])
    (let* ( (segments (tm-seg-register tm))
            (sorted (lambda (counts flatxs threshold)
                      (list-sort (lambda (sega segb)
                                   (if (fx=? (seg-cellx sega) (seg-cellx segb))
                                       (fx<? (seg-flatx sega) (seg-flatx segb))
                                       (fx<? (seg-cellx sega) (seg-cellx segb))))
                        (do ((i 0 (add1 i)) 
                             (segs '() (if (fx>=? (vector-ref counts i) threshold)
                                         (cons (vector-ref segments (vector-ref flatxs i)) segs)
                                         segs)))
                            ((fx=? i (vector-length counts)) segs))))))
      (tm-active-segments-set!   tm                                            ;; 67. activeSegments(t).add(segment)
          (sorted active-counts active-flatxs (tm-activation-threshold tm)))   ;; 66. if numActiveConnected ≥ ACTIVATION_THRESHOLD then
      (tm-matching-segments-set! tm                                            ;; 70. matchingSegments(t).add(segment)
          (sorted matching-counts matching-flatxs (tm-min-threshold tm)))      ;; 69. if numActivePotential ≥ LEARNING_THRESHOLD then
      (tm-num-active-pot-syns-for-seg-set! tm num-active-potential)))          ;; 72. numActivePotentialSynapses(t, segment) = numActivePotential
  (when learn
    (let loop ((actsegs (tm-active-segments tm)))
      (cond [ (null? actsegs) (tm-iteration-set! tm (add1 (tm-iteration tm))) ]
            [ else (seg-last-used-set! (car actsegs) (tm-iteration tm))
                   (loop (cdr actsegs)) ]))))
                                                                                            ;
;; --- Supporting functions ---
                                                                                            ;
(define (least-used-cell tm cells)       ;; TM (listof CellX) -> CellX
  ;; produce a cell with fewest segments                                       ;; 73. function leastUsedCell(column)
  (let loop ((cells cells) (candidates '()) (fewest-segs (greatest-fixnum)))   ;; 75/79,78. for cell in column.cells; leastUsedCells = []
    (if (null? cells) (list-ref candidates (random (length candidates)))       ;; 83. return chooseRandom(leastUsedCells)
        (let ((n-segs (length (vector-ref (tm-cells tm) (car cells)))))
          (cond [ (fx<? n-segs fewest-segs)                                    ;; 76. fewestSegments = min(fewestSegments, cell.segments.length)
                  (loop (cdr cells) (list (car cells)) n-segs) ]
                [ (fx=? n-segs fewest-segs)                                    ;; 80. if cell.segments.length == fewestSegments then
                  (loop (cdr cells) (cons (car cells) candidates) n-segs) ]    ;; 81.   leastUsedCells.add(cell)
                [ else (loop (cdr cells) candidates fewest-segs) ] )))))

  #;( [vector-set! (tm-cells TM22) 0 '(SEG00)]
      [vector-set! (tm-cells TM22) 1 '()]
      [expect ([least-used-cell TM22 '(0 1)] 1 )] )
                                                                                            ;
(define (best-matching-segment tm column ;; TM ColX {Seg} -> {Seg} {Seg}
          column-matching-segments)                                            ;; 84. function bestMatchingSegment(column)
  ;; step through matching segments for this column; return next and seg with most synapses 
  (let loop ( (segments column-matching-segments)                              ;; 87. for segment in segmentsForColumn(column, matchingSegments(t-1))
              (bestseg  column-matching-segments)                              ;; 85. bestMatchingSegment = None
              (best-score -1))                                                 ;; 86. bestScore = -1
    (cond [ (fx=? column (segs->colx tm segments))
            (let ((naps (hashtable-ref (tm-num-active-pot-syns-for-seg tm)
                                       (seg-flatx (car segments)) 0)))
              (if (fx>? naps best-score)                                       ;; 88. if numActivePotentialSynapses(t-1, segment) > bestScore then
                  (loop (cdr segments) segments naps)
                  (loop (cdr segments) bestseg best-score))) ]                 ;; 89-90. bestMatchingSegment = segment; bestScore = nAPS(t-1, segment)
          [ else (values segments bestseg) ])))                                ;; 92. return bestMatchingSegment
                                                                                            ;
(define (grow-synapses tm                ;; TM Seg Nat (listof CellX) ->
          segment n-desired-new-synapses prev-winner-cells)                    ;; 93. function growSynapses(segment, newSynapseCount)
  ;; create synapses from winners on segment, replacing low-permanence ones as needed
  (let* ( (synapses     (seg-synapses segment))
          (num-synapses (synapses-length synapses))
          (candidates   (let loop ((cs '()) (pwc prev-winner-cells))           ;; 94. candidates = copy(winnerCells(t-1))
                          (cond [ (null? pwc) cs ]
                                [ (memv-prex? (car pwc) synapses)              ;; 99-104. [omit from candidates if alreadyConnected]
                                       (loop cs (cdr pwc)) ]
                                [ else (loop (cons (car pwc) cs) (cdr pwc)) ])))
          (n-actual (fxmin n-desired-new-synapses (length candidates)))
          (overrun  (fx- (fx+ num-synapses n-actual) (tm-max-synapses-per-segment tm))))
    (when (fxpositive? overrun)
      (destroy-min-permanence-synapses tm segment overrun prev-winner-cells))
    (let* ( (num-synapses (synapses-length synapses))
            (n-actual (fxmin n-actual (fx- (tm-max-synapses-per-segment tm) num-synapses)))
            (new-prexs (vector-sample (list->vector candidates) n-actual)))    ;; 96. presynapticCell = chooseRandom(candidates)
      (seg-synapses-set! segment
        (build-synapses (fx+ num-synapses n-actual)                            ;; 95. while candidates.length > 0 and newSynapseCount > 0
          (lambda (sx)
            (if (fx<? sx num-synapses)
                (synapses-ref synapses sx)            ;; copy over retained synapses
                (let ((new-prex (vector-ref new-prexs (fx- sx num-synapses))))
                  (add-to-pre-index tm new-prex segment)
                  (make-synapse new-prex (tm-initial-permanence tm))))))))))   ;; 105. newSynapse = createNewSynapse(segment, presynapticCell, INITIAL_PERMANENCE)

  #;( [grow-synapses TM22 SEG00 2 '(2 3)]
      (let ((ms (lambda (prex) [make-synapse prex [tm-initial-permanence TM22]])))
        [expect ([or (equal? [seg-synapses SEG00] (synapses: (ms 2) (ms 3)))
                     (equal? [seg-synapses SEG00] (synapses: (ms 3) (ms 2)))] #t )]) )
                                                                                            ;
(define (skip-col tm column segments)    ;; TM ColX (listof Seg) -> (listof Seg)
  ;; step to next column's segments
  (let loop ((segments segments))
    (if (fx=? column (segs->colx tm segments))
        (loop (cdr segments))
        segments)))
                                                                                            ;
;; --- Exported functions ---
                                                                                            ;
(define (tm-compute tm                   ;; TM (listof ColumnX) Boolean ->
          active-columns learn)
  ;; one iteration
  (activate-cells tm (list-sort < active-columns) learn)
  (activate-dendrites tm learn))
                                                                                            ;
(define (tm-reset tm)                    ;; TM ->
  (tm-active-cells-set!      tm '())
  (tm-winner-cells-set!      tm '())
  (tm-active-segments-set!   tm '())
  (tm-matching-segments-set! tm '()))
                                                                                            ;
(define (tm-get-active-cols tm)          ;; TM -> (listof ColX)
  (let loop ((previous-colx -1) (active-cols '()) (cells (tm-active-cells tm)))
    (cond [ (null? cells) active-cols ]
          [ (fx=? (cellx->colx tm (car cells)) previous-colx) 
              (loop previous-colx active-cols (cdr cells)) ]
          [ else (loop (cellx->colx tm (car cells)) 
                       (cons (cellx->colx tm (car cells)) active-cols) (cdr cells)) ])))
                                                                                                ;
(define (tm-get-predictive-cols tm)      ;; TM -> (listof ColX)
  (let loop ((previous-colx -1) (predictive-cols '()) (segments (tm-active-segments tm)))
    (let ((this-colx (segs->colx tm segments)))
      (cond [ (null? segments) predictive-cols ]
            [ (fx=? this-colx previous-colx) 
                (loop previous-colx predictive-cols (cdr segments)) ]
            [ else (loop this-colx (cons this-colx predictive-cols) (cdr segments)) ]))))
    
)