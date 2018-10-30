#!r6rs

;; ========= HTM-scheme Temporal Memory Copyright 2017 Roger Turner. =========
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Based on code from Numenta Platform for Intelligent Computing (NuPIC) ;;
  ;; which is Copyright (C) 2014-2016, Numenta, Inc.                       ;;
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

  ;; Translated from numenta nupic/.../temporal_memory.py and connections.py -- 
  ;; see comments there for descriptions of functions and parameters.
  ;; Numbered comments in the "Temporal Memory Algorithm" section echo, where
  ;; possible, pseudocode in Numenta BaMI Temporal-Memory-Algorithm-Details.pdf 
  ;; revision 0.5
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.

(library (HTM-scheme HTM-scheme algorithms temporal_memory)
                                                                                            ;
(export
  tm-num-columns
  tm-cells-per-column
  tm-initial-permanence
  tm-max-new-synapse-count
  tm-permanence-increment
  tm-permanence-decrement
  tm-predicted-segment-decrement
  tm-next-flatx
  tm-active-cells
  tm-active-cells-set!
  tm-winner-cells
  tm-winner-cells-set!
  tm-active-segments
  tm-active-segments-set!
  tm-matching-segments
  tm-matching-segments-set!
  tm-max-synapses-per-segment
  tm-num-cells
  tm-cells
  tm-pre-index
  tm-iteration
  tm-iteration-set!
  tm-num-active-pot-syns-for-seg
  punish-predicted-column
  compute-activity
  sorted-segs
  adapt-segment
  create-segment
  least-used-cell
  skip-col
  cellx->colx
  segs->cellx
  segs->colx
  add-to-pre-index
  colxs->colx
  seg-cellx
  seg-flatx
  seg-synapses
  seg-synapses-set!
  seg-last-used-set!
  syn-perm
  synapses-ref
  build-synapses
  in-synapses?
  syn-prex
  synapses-length
  clip-max
  clip-min
  (rename
    (tm                  tm:tm)
    (make-tm             tm:construct)
    (perm<-              tm:permanence)
    (make-synapse        tm:synapse)
    (make-seg            tm:segment)
    (compute             tm:compute)
    (reset               tm:reset)
    (get-active-cols     tm:get-active-cols)
    (get-predictive-cols tm:get-predictive-cols)))
                                                                                            ;    
(import 
  (rnrs)
  (HTM-scheme HTM-scheme algorithms htm_prelude))

;; === Temporal Memory Types ===
                                                                                            ;
;; Boolean, Number, Integer, Fixnum, (listof X), (vectorof X) ... = Scheme types
;; X, Y, Z      = type parameters (arbitrary types in function type specification etc)
;; X Y -> Z     = function with argument types X, Y and result type Z
;; {X}          = abbreviation for (listof X)
;; Nat          = natural number (including zero) (Scheme Fixnum or exact Integer)
;; Fix10k       = Fixnum interpreted as number with 4 decimal places, ie 10000 = 1.0
;; Permanence   = Fix10k [0 - MAXPERM], interpreted as 0.0 - 0.9999
;; Synapse      = Fixnum, interpreted as CellX * 10000 + Permanence 
;; Synapses     = (vectorof Synapse) [optionally fxvector in Chez Scheme]
;; Segment      = Record: CellX, FlatX, and Synapses
;; FlatX        = Nat, segment sequence number
;; Cell         = (listof Segment)
;; CellX        = Nat [0 - (div (expt 2 (fixnum-width)) 10000)] cell index
;; CellVecOf    = Vector indexed by CellX
;; ColX         = Nat [0 - MAXCOL], column index of cell
;; TM           = Record: tm parameters, cells

;; === Parameters and Data ===
                                                                                            ;
(define %max-perm% 9999)                 ;; Fixnum
(define %min-perm% 0)                    ;; Fixnum
                                                                                            ;
(define (clip-max perm)                  ;; Permanence -> Permanence
  (fxmin %max-perm% perm))
                                                                                            ;
(define (clip-min perm)                  ;; Permanence -> Permanence
  (fxmax %min-perm% perm))
                                                                                            ;
(define (perm<- x)                       ;; Number[0.0-.9999] -> Permanence
  (clip-max (clip-min (int<- (* x x10k)))))
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
    (mutable next-flatx)                 ;; Nat: next available index in seg-index
    (mutable free-flatx)                 ;; (listof Nat): indices available for re-use
    (mutable active-cells)               ;; (vectorof CellX)
    (mutable winner-cells)               ;; (listof CellX)
    (mutable active-segments)            ;; (listof Segment)
    (mutable matching-segments)          ;; (listof Segment)
    (mutable iteration)                  ;; Nat
    prediction-fail-boost                ;; Nat: predicted-segment-decrement boost range
    num-columns                          ;; Nat
    num-cells                            ;; Nat
    cells                                ;; (CellVecOf {Segment})
    pre-index                            ;; (CellVecOf {FlatX})
    (mutable seg-index)                  ;; (vectorof Segment)
    (mutable num-active-pot-syns-for-seg);; (vectorof Nat)
    )
  (protocol
    (lambda (new)
      (lambda (cd cpc . kwargs)            ;; (listof Nat) Nat (listof KWarg) -> TM
        ;; produce temporal memory instance with cd columns, cpc cells-per-column, and
        ;; defaults for parameters not specified by args; seg-index receives initial
        ;; allocation, extended in create-segment as needed
        (let* ( (num-columns (apply * cd))
                (num-cells   (* num-columns cpc)))
          (apply new (append 
                       (list cd cpc) 
                       (key-word-args
                          (append 
                            `([num-columns                 . ,num-columns]
                              [num-cells                   . ,num-cells]
                              [cells                       . ,(make-vector num-cells '())]
                              [pre-index                   . ,(make-vector num-cells '())]
                              [seg-index                   . ,(make-vector num-columns)])
                            kwargs)
                          tm-defaults))))))))
                                                                                            ;
(define tm-defaults                      ;; (listof KWarg)
  `(
    [activation-threshold        . 13]
    [initial-permanence          . ,(perm<- 0.21)]
    [connected-permanence        . ,(perm<- 0.5)]
    [min-threshold               . 10]
    [max-new-synapse-count       . 20]
    [permanence-increment        . ,(perm<- 0.1)]
    [permanence-decrement        . ,(perm<- 0.1)]
    [predicted-segment-decrement . ,(perm<- 0.004)]
    [max-segments-per-cell       . 255]
    [max-synapses-per-segment    . 255]
    [next-flatx                  . 0]
    [free-flatx                  . ()]
    [active-cells                . #()]
    [winner-cells                . ()]
    [active-segments             . ()]
    [matching-segments           . ()]
    [iteration                   . 0]
    [prediction-fail-boost       . 1]
    [num-columns                 . 0]
    [num-cells                   . 0]
    [cells                       . #()]
    [pre-index                   . #()]
    [seg-index                   . #()]
    [num-active-pot-syns-for-seg . #()]
    ))
                                                                                            ;
(define-record-type seg                  ;; Segment
  (fields
    cellx                                ;; Nat: index of cell that this is a segment of
    flatx                                ;; Nat: index of this in the register of segments
    (mutable last-used)                  ;; Nat: iteration number last time segment active
    (mutable last-punished)              ;; Nat: iteration number when last punished
    (mutable synapses)                   ;; (vectorof Synapse)
    ))

;; === Temporal Memory Algorithm ===
                                                                                            ;
;; --- Main functions: see NuPIC temporal_memory.py ---
                                                                                            ;
(define (compute tm active-columns learn);; TM (listof ColX) Boolean ->
  ;; Perform one time step of the Temporal Memory algorithm.
  (activate-cells tm (list-sort fx<? active-columns) learn)
  (activate-dendrites tm learn))
                                                                                            ;
(define (activate-cells tm               ;; TM (listof ColX) Boolean ->
          active-columns learn)
  ;; Calculate the active cells, using the current active columns and dendrite
  ;; segments. Grow and reinforce synapses.
  ;; Step through columns in merge of actcols/actsegs/matsegs lists, building active/winner lists; 
  ;; the activate/burst/punish functions process one column and return list for next column.
  ;; This blends the Python groupby2 logic into these functions; actcells/wincells are current,
  ;; tm-active/winner-cells are the values from the previous time step.
  ;; Numbered comments echo pseudocode in BaMI Temporal-Memory-Algorithm-Details.pdf revision 0.5
  (let loop ( (actcols active-columns)                                         ;; 1. for column in columns
              (actsegs (tm-active-segments tm)) 
              (matsegs (tm-matching-segments tm)) 
              (actcells '()) 
              (wincells '()))
    (let ((actcols-colx (colxs->colx tm actcols))
          (actsegs-colx (segs->colx  tm actsegs))
          (matsegs-colx (segs->colx  tm matsegs)))
      (let ((column (fxmin actcols-colx actsegs-colx matsegs-colx)))
        (if (fx<? column (tm-num-columns tm))
          (if (fx=? column actcols-colx)                                       ;; 2. if column in activeColumns(t) then
              (if (fx=? column actsegs-colx)                                   ;; 3.   if count(segmentsForColumn(column, activeSegments(t-1))) > 0 then
                  (let-values ([(nextactseg actcells wincells)                 ;; 4.     activatePredictedColumn(column)
                      (activate-predicted-column tm column actsegs actcells wincells learn)])
                    (loop actcols nextactseg matsegs actcells wincells))
                  (let-values ([(nextmatseg actcells wincells)                 ;; 6.   else burstColumn(column)
                      (burst-column tm column matsegs actcells wincells learn)])
                    (loop (if (null? actcols) actcols (cdr actcols))
                          actsegs nextmatseg actcells wincells)))
              (if (fx=? column matsegs-colx)                                   ;; 8.  else if count(segmentsForColumn(column, matchingSegments(t-1))) > 0 then
                  (let ((nextmatseg (if learn                                  ;; 50.   if LEARNING_ENABLED
                                      (punish-predicted-column                 ;; 9.      punishPredictedColumn(column)
                                        tm column matsegs (tm-active-cells tm)
                                        (tm-predicted-segment-decrement tm))
                                      (skip-col tm column matsegs))))
                    (loop actcols actsegs nextmatseg actcells wincells))
                  (loop actcols (skip-col tm column actsegs) matsegs actcells wincells)))
          ;; all active cols and segs handled: set prev active/winner cells for next iteration
          (begin
            (tm-active-cells-set! tm (list->vector actcells))
            (vector-sort! < (tm-active-cells tm))
            (tm-winner-cells-set! tm wincells)))))))
                                                                                            ;
(define (activate-dendrites tm learn)    ;; TM Boolean ->
  ;; Calculate dendrite segment activity, using the current active cells.
  (let-values ([(actsegs potsegs num-active-potential) 
                (compute-activity 
                  tm (tm-active-cells tm) (tm-pre-index tm) 0 '())])           ;; 66. if numActiveConnected ≥ ACTIVATION_THRESHOLD then
    (tm-active-segments-set! tm (sorted-segs actsegs))                         ;; 67. activeSegments(t).add(segment)
                                                                               ;; 69. if numActivePotential ≥ LEARNING_THRESHOLD then
    (tm-matching-segments-set! tm (sorted-segs potsegs))                       ;; 70. matchingSegments(t).add(segment)
    (tm-num-active-pot-syns-for-seg-set! tm num-active-potential))             ;; 72. numActivePotentialSynapses(t, segment) = numActivePotential
  (when learn
    (for-each
      (lambda (segment)
        (seg-last-used-set! segment (tm-iteration tm)))
      (tm-active-segments tm))
    (tm-iteration-set! tm (add1 (tm-iteration tm)))))
                                                                                            ;
(define (reset tm)                       ;; TM ->
  ;; Indicates the start of a new sequence. Clears any predictions and makes sure
  ;; synapses don't grow to the currently active cells in the next time step.
  (tm-active-cells-set!      tm '#())
  (tm-winner-cells-set!      tm '())
  (tm-active-segments-set!   tm '())
  (tm-matching-segments-set! tm '()))
                                                                                            ;
(define (activate-predicted-column tm    ;; TM ColX {Segment} {CellX} {CellX} Boolean -> {Segment} {CellX} {CellX}
          column active-segs actcells wincells learn)                          ;; 10. function activatePredictedColumn(column)
  ;; Determines which cells in a predicted column should be added to winner cells
  ;; list, and learns on the segments that correctly predicted this column.
  (let loop ( (segments  active-segs) 
              (actcells  actcells) 
              (wincells  wincells) 
              (prev-cell (tm-num-cells tm)))
    (if (fx=? (segs->colx tm segments) column)                                 ;; 11. for segment in segmentsForColumn(column, activeSegments(t-1))
        (let ((segment (car segments)))
          (when learn                                                          ;; 15. if LEARNING_ENABLED
            (learn-on-cell tm segment))
          (let ((this-cell (seg-cellx segment)))
            (if (fx=? this-cell prev-cell)
                (loop (cdr segments) actcells wincells prev-cell)
                (loop (cdr segments) 
                      (cons this-cell actcells)                                ;; 12. activeCells(t).add(segment.cell)
                      (cons this-cell wincells)                                ;; 13. winnerCells(t).add(segment.cell)
                      this-cell))))
        ;; else: end of segments for this column, return remaining segments and updated cell lists
        (values segments actcells wincells))))
                                                                                            ;
(define (burst-column tm column          ;; TM ColX {Segment} {CellX} {CellX} Boolean -> {Segment} {CellX} {CellX}
          column-matching-segments actcells wincells learn)                    ;; 25. function burstColumn(column)
  ;; Activates all of the cells in an unpredicted active column, chooses a winner
  ;; cell, and, if learning is turned on, learns on one segment, growing a new
  ;; segment if necessary.
  (let* ( (start (fx* column (tm-cells-per-column tm)))
          (cells-for-column (build-list (tm-cells-per-column tm)               ;; 26. for cell in column.cells
                                        (lambda (x) (fx+ x start)))))
    (if (fx=? column (segs->colx tm column-matching-segments))                 ;; 29. if segmentsForColumn(column, matchingSegments(t-1)).length > 0
        (let-values ([(nextmatseg best-match)                                  ;; 30.   learningSegment = bestMatchingSegment(column)
                      (best-matching-segment tm column column-matching-segments)])
          (when learn                                                          ;; 39. if LEARNING_ENABLED
            (learn-on-cell tm (car best-match)))                               ;; 40-48. [use learn-on-cell]
          (values nextmatseg 
                  (append cells-for-column actcells)                           ;; 27. activeCells(t).add(cell)
                  (cons (seg-cellx (car best-match)) wincells)))               ;; 37/31. winnerCells(t).add(learningSegment.cell)
        ;; else: no matching segments for column
        (let ((winner-cell (least-used-cell tm cells-for-column)))             ;; 33. else winnerCell = leastUsedCell(column)
          (when learn                                                          ;; 34.   if LEARNING_ENABLED
            (let ((n-grow-exact (fxmin (tm-max-new-synapse-count tm)
                                       (length (tm-winner-cells tm)))))
              (when (fxpositive? n-grow-exact)                                 ;; 35.   learningSegment = growNewSegment(winnerCell)
                (let ((lseg (create-segment tm winner-cell)))                  ;; 48.   growSynapses(learningSegment, newSynapseCount)
                  (grow-synapses tm lseg n-grow-exact 
                                 (tm-winner-cells tm) (tm-pre-index tm))))))
          (values column-matching-segments 
                  (append cells-for-column actcells)                           ;; 27. activeCells(t).add(cell)
                  (cons winner-cell wincells))))))                             ;; 37. winnerCells(t).add(winnerCell)
                                                                                            ;
(define (punish-predicted-column         ;; TM ColX {Segment} (vectorof CellX) -> {Segment}
          tm colx column-matching-segments prev-active-cells decrement)        ;; 49. function punishPredictedColumn(column)
  ;; Punishes the Segments that incorrectly predicted a column to be active.
  (if (null? column-matching-segments) column-matching-segments
    (let loop ((segments column-matching-segments))                            ;; 51. for segment in segmentsForColumn(column, matchingSegments(t-1))
      (if (fx=? colx (segs->colx tm segments))
        (let ((segment   (car segments)))
          (when (fxpositive? decrement)                                        ;; 52-54. [use adapt-segment to decrement permanences]
            (adapt-segment tm segment prev-active-cells (tm-cells tm)
                ;; scale up predicted-segment-decrement on repeated misprediction
                ;; within prediction-fail-boost iterations (*not in NuPIC*)
                (lambda (p)
                  (let ((since-punished (fx- (tm-iteration tm) (seg-last-punished segment))))
                    (clip-min (fx- p (if (fx<? since-punished (tm-prediction-fail-boost tm))
                                         (fx* (fx- (add1 (tm-prediction-fail-boost tm))
                                                   since-punished) 
                                              decrement)
                                         decrement)))))
                0))
          (seg-last-punished-set! segment (tm-iteration tm))
          (loop (cdr segments)))
        ;; else: end of segments for this column, return remaining segments
        segments))))
                                                                                            ;
(define  create-segment                  ;; TM CellX [(CellVecOf {Segment})] -> Seg
  (case-lambda 
  [(tm cellx)
    ;; default connections
    (create-segment tm cellx (tm-cells tm))]
  [(tm cellx connections)
    ;; produce a new segment on the cell, updating the index of segments
    (let ((segments (vector-ref connections cellx)))
      (when (fx>=? (length segments) (tm-max-segments-per-cell tm))
        (let loop ((segs segments) (oldest-seg #f) (last-used (greatest-fixnum)))
          (cond [ (null? segs)
                  (vector-set! connections cellx (remq oldest-seg segments))
                  (destroy-segment tm (car oldest-seg)) ]
                [ (fx<? (seg-last-used (car segs)) last-used)
                  (loop (cdr segs) segs (seg-last-used (car segs))) ]
                [ else (loop (cdr segs) oldest-seg last-used) ])))
      (let* ( (new-flatx (if (null? (tm-free-flatx tm))
                            (let ((next (tm-next-flatx tm)))
                              (when (fx=? next (vector-length (tm-seg-index tm)))
                                (tm-seg-index-set! tm (vector-extend (tm-seg-index tm))))
                              next)  
                            (car (tm-free-flatx tm))))
              (segment (make-seg cellx new-flatx (tm-iteration tm)
                                 (fx- (tm-prediction-fail-boost tm)) (make-synapses 0))))
        (vector-set! (tm-seg-index tm) new-flatx segment)
        (if (null? (tm-free-flatx tm))
          (tm-next-flatx-set! tm (add1 new-flatx))
          (tm-free-flatx-set! tm (cdr (tm-free-flatx tm))))
        (vector-set! connections cellx (cons segment segments))
        segment))]))
                                                                                            ;
(define (destroy-min-permanence-synapses ;; TM Segment Nat (listof CellX) -> 
          tm segment n-destroy exclude-cells)
  ;; drop n-destroy lowest-permanence synapses, but retaining any on exclude list
  (let-values ([(keepers destroy-candidates)
                (partition (lambda (x) (memv (syn-prex x) exclude-cells))
                           (synapses->list (seg-synapses segment)))])
    (let* ( (destroy-candidates 
              (list-sort 
                (lambda (x y) (fx<? (syn-perm x) (syn-perm y))) 
                destroy-candidates))
            (survivors (list-tail destroy-candidates (fxmin n-destroy (length destroy-candidates)))))
      (seg-synapses-set! segment (list->synapses (append survivors keepers))))))
                                                                                            ;
(define (least-used-cell tm cells)       ;; TM (listof CellX) -> CellX
  ;; produce a cell with fewest segments; break ties randomly                  ;; 73. function leastUsedCell(column)
  (let loop ((cells cells) (candidates '()) (fewest-segs (greatest-fixnum)))   ;; 75/79,78. for cell in column.cells; leastUsedCells = []
    (if (null? cells) 
        (list-ref candidates (random (length candidates)))                     ;; 83. return chooseRandom(leastUsedCells)
        (let ((n-segs (length (vector-ref (tm-cells tm) (car cells)))))
          (cond [ (fx<? n-segs fewest-segs)                                    ;; 76. fewestSegments = min(fewestSegments, cell.segments.length)
                  (loop (cdr cells) (list (car cells)) n-segs) ]
                [ (fx=? n-segs fewest-segs)                                    ;; 80. if cell.segments.length == fewestSegments then
                  (loop (cdr cells) (cons (car cells) candidates) n-segs) ]    ;; 81.   leastUsedCells.add(cell)
                [ else (loop (cdr cells) candidates fewest-segs) ] )))))
                                                                                            ;
(define (grow-synapses                   ;; TM Segment Nat (listof CellX) [(CellVecOf {FlatX})] ->
          tm segment n-desired-new-synapses prev-winner-cells pre-index)       ;; 93. function growSynapses(segment, newSynapseCount)
  ;; create synapses from winners on segment, replacing low-permanence ones as needed
  (let* ( (synapses     (seg-synapses segment))
          (num-synapses (synapses-length synapses))
          (candidates   (let loop ((cs '()) (pwc prev-winner-cells))           ;; 94. candidates = copy(winnerCells(t-1))
                          (cond [ (null? pwc) cs ]
                                [ (in-synapses? (car pwc) synapses)            ;; 99-104. [omit from candidates if alreadyConnected]
                                    (loop cs (cdr pwc)) ]
                                [ else (loop (cons (car pwc) cs) (cdr pwc)) ])))
          (n-actual     (fxmin n-desired-new-synapses (length candidates)))
          (max-synapses-per-segment 
                        (if (positive? (tm-max-synapses-per-segment tm))
                          (tm-max-synapses-per-segment tm)
                          (fxdiv (greatest-fixnum) 2)))
          (overrun      (fx- (fx+ num-synapses n-actual) max-synapses-per-segment)))
    (when (fxpositive? overrun)
      (destroy-min-permanence-synapses tm segment overrun prev-winner-cells))
    (let* ( (num-synapses (synapses-length synapses))
            (n-actual (fxmin n-actual (fx- max-synapses-per-segment num-synapses))))
      (when (positive? n-actual)
        (let ((new-prexs (vector-sample (list->vector candidates) n-actual)))  ;; 96. presynapticCell = chooseRandom(candidates)
          (seg-synapses-set! segment
            (build-synapses (fx+ num-synapses n-actual)                        ;; 95. while candidates.length > 0 and newSynapseCount > 0
              (lambda (sx)
                ;; copy retained synapses, make n-actual new ones
                (if (fx<? sx num-synapses)
                    (synapses-ref synapses sx)
                    (let ((new-prex (vector-ref new-prexs (fx- sx num-synapses))))
                      (add-to-pre-index new-prex segment pre-index)
                      (make-synapse new-prex (tm-initial-permanence tm)))))))  ;; 105. newSynapse = createNewSynapse(segment, presynapticCell, INITIAL_PERMANENCE)
          ;; keep synapses sorted for binary search in compute-activity
          (vector-sort! fx<? (seg-synapses segment)))))))
                                                                                            ;
(define  adapt-segment                   ;; TM {Segment} (vectorof CellX) [(CellVecOf {Segment}) [(Perm -> Perm) Permanence]] ->
  ;; Updates synapses on segment. Strengthens active synapses; weakens inactive synapses.
  ;; Remove synapse on zero permanence, destroy segment if no synapses left.
  (case-lambda
  [(tm segment active-input)
    ;; default connections
    (adapt-segment tm segment active-input (tm-cells tm))]
  [(tm segment active-input connections)
    ;; use increment/decrement parameter values
    (adapt-segment tm segment active-input connections
                   (lambda (p) (clip-max (fx+ p (tm-permanence-increment tm))))
                   (tm-permanence-decrement tm))]
  [(tm segment active-input connections increment-proc permanence-decrement)
    ;; general case: increment supplied by proc
    (let ((synapses    (seg-synapses segment))
          (input-right (fx- (vector-length active-input) 1)))
      (let build-s-t-d ((sx (fx- (synapses-length synapses) 1))                ;; 16. for synapse in segment.synapses
                        (synapses-to-destroy '()))
        (if (negative? sx)
            (cond [ (null? synapses-to-destroy) ]
                  [ (fx=? (length synapses-to-destroy) (synapses-length synapses))
                      (vector-set! connections (seg-cellx segment)
                        (remq segment (vector-ref connections (seg-cellx segment))))
                      (destroy-segment tm segment) ]
                  [ else (seg-synapses-set! segment 
                           (prune-synapses synapses synapses-to-destroy)) ])
            (let* ( (synapse (synapses-ref synapses sx))
                    (prex    (syn-prex synapse))
                    (permanence
                      ;; binary search for index of presynaptic cell in active-input
                      (let search ((left 0) (right input-right))
                        (if (fx>? left right)                                  ;; 20. else synapse.permanence -= PERMANENCE_DECREMENT
                          (clip-min (fx- (syn-perm synapse) permanence-decrement))
                          (let ((mid (fxdiv (fx+ left right) 2)))
                            (cond                                              ;; 17. if synapse.presynapticCell in activeCells(t-1) then
                              [ (fx<? prex (vector-ref active-input mid)) 
                                  (search left (fx- mid 1)) ]
                              [ (fx<? (vector-ref active-input mid) prex) 
                                  (search (add1 mid) right) ]
                              [ else (increment-proc (syn-perm synapse))])))))) ;; 18.   synapse.permanence += PERMANENCE_INCREMENT
              (if (zero? permanence)
                  ;; build synapses-to-destroy indices as sorted list
                  (build-s-t-d (fx- sx 1) (cons sx synapses-to-destroy))
                  (begin
                    (synapses-set! synapses sx (make-synapse prex permanence))
                    (build-s-t-d (fx- sx 1) synapses-to-destroy)))))))]))
                                                                                            ;
(define (learn-on-cell tm segment)       ;; TM Segment ->
  ;; adjust permanences and add new synapses to segment
  (adapt-segment tm segment (tm-active-cells tm))
  (let ((n-grow-desired (fx- (tm-max-new-synapse-count tm)                     ;; 22. newSynapseCount = (SYNAPSE_SAMPLE_SIZE -
                             (vector-ref (tm-num-active-pot-syns-for-seg tm)   ;; 23.   numActivePotentialSynapses(t-1, segment))
                                         (seg-flatx segment)))))
    (when (fxpositive? n-grow-desired)                                         ;; 24. growSynapses(segment, newSynapseCount)
      (grow-synapses tm segment n-grow-desired 
                     (tm-winner-cells tm) (tm-pre-index tm)))))
                                                                                            ;
(define (get-active-cols tm)             ;; TM -> (listof ColX)
  (let ((cells (tm-active-cells tm)))
    (let loop ((previous-colx -1) (active-cols '()) (cx 0))
      (if (fx<? cx (vector-length cells))
        (let ((colx (cellx->colx tm (vector-ref cells cx))))
          (if (fx=? colx previous-colx) 
              (loop previous-colx active-cols (add1 cx))
              (loop colx (cons colx active-cols) (add1 cx))))
        active-cols))))
                                                                                            ;
(define (get-predictive-cols tm)         ;; TM -> (listof ColX)
  (let loop ((previous-colx -1) (predictive-cols '()) (segments (tm-active-segments tm)))
    (let ((this-colx (segs->colx tm segments)))
      (cond [ (null? segments) (reverse predictive-cols) ]
            [ (fx=? this-colx previous-colx) 
                (loop previous-colx predictive-cols (cdr segments)) ]
            [ else (loop this-colx (cons this-colx predictive-cols) (cdr segments)) ]))))

;; === Supporting Functions ===
                                                                                            ;
;; --- Connections: see NuPIC connections.py ---
                                                                                            ;
(define (compute-activity                ;; TM (vectorof CellX) (CellVecOf {FlatX}) -> {Segment} {Segment} (vectorof Nat)
          tm active-presynaptic-cells pre-index reduced-threshold reduced-threshold-cells)
  ;; produce segments with connected/potential synapses, and counts, for active cells
  (let* (
    (napsfs    (make-vector (tm-next-flatx tm) 0))    ;; "num-active-potential-synapses-for-segment"
    (nacsfs    (make-vector (tm-next-flatx tm) 0))    ;; "num-active-connected-synapses-for-segment"
    (threshold (tm-connected-permanence tm))
    (segments  (tm-seg-index tm))
    (actsegs   '())
    (potsegs   '())
    (actthresh (tm-activation-threshold tm))
    (potthresh (tm-min-threshold tm)))
    (define (update-actsegs segx synapse)
      (when (fx>=? (syn-perm synapse) threshold)                               ;; 60. if synapse.permanence ≥ CONNECTED_PERMANENCE then
        (let ((segment    (vector-ref segments segx))
              (new-nacsfs (add1 (vector-ref nacsfs segx))))
          (if (null? reduced-threshold-cells)
              (when (fx=? new-nacsfs actthresh)                                ;; 61.   numActiveConnected += 1
                (unless (memq segment actsegs)
                  (set! actsegs (cons segment actsegs))))
              (when (or (and (fx=? new-nacsfs reduced-threshold)
                             (memv (seg-cellx segment) reduced-threshold-cells))
                        (fx=? new-nacsfs actthresh))
                (unless (memq segment actsegs)
                  (set! actsegs (cons segment actsegs)))))
          (vector-set! nacsfs segx new-nacsfs))))
    (define (update-potsegs segx syn-low syn-high)
      (let ((synapses (seg-synapses (vector-ref segments segx))))
        (let search ((left 0) (right (fx- (synapses-length synapses) 1)))
          (unless (fx>? left right)
            (let* ( (mid (fxdiv (fx+ left right) 2))
                    (synapse (synapses-ref synapses mid))) 
              (cond 
                [ (fx<? synapse  syn-low) (search (add1 mid) right) ]
                [ (fx<? syn-high synapse) (search left (fx- mid 1)) ]
                [ else                                                         ;; 63. if synapse.permanence ≥ 0 then
                  (let ((new-napsfs (add1 (vector-ref napsfs segx))))
                    (when (fx=? new-napsfs potthresh)
                      (unless (memq (vector-ref segments segx) potsegs)
                        (set! potsegs (cons (vector-ref segments segx) potsegs))))
                    (vector-set! napsfs segx new-napsfs))                      ;; 64.   numActivePotential += 1
                  (update-actsegs segx synapse) ]))))))
    (vector-for-each                                                           ;; 59. if synapse.presynapticCell in activeCells(t) then
      (lambda (cellx)
        (let* ((syn-low  (make-synapse cellx %min-perm%))
               (syn-high (fx+ syn-low %max-perm%)))
          (for-each                                                            ;; 55. for segment in segments
            (lambda (segx)
              (update-potsegs segx syn-low syn-high))
            (vector-ref pre-index cellx))))
      active-presynaptic-cells)
    (values actsegs potsegs napsfs)))
                                                                                            ;
;; --- Segments ---
                                                                                            ;
(define (best-matching-segment tm column ;; TM ColX {Segment} -> {Segment} {Segment}
          column-matching-segments)                                            ;; 84. function bestMatchingSegment(column)
  ;; step thru matching segs for this col; return next/seg with most synapses 
  (let loop ( (segments column-matching-segments)                              ;; 87. for segment in segmentsForColumn(column, matchingSegments(t-1))
              (bestseg  column-matching-segments)                              ;; 85. bestMatchingSegment = None
              (best-score -1))                                                 ;; 86. bestScore = -1
    (if (fx=? column (segs->colx tm segments))
        (let ((naps (vector-ref (tm-num-active-pot-syns-for-seg tm)
                                       (seg-flatx (car segments)))))
          (if (fx>? naps best-score)                                           ;; 88. if numActivePotentialSynapses(t-1, segment) > bestScore then
              (loop (cdr segments) segments naps)
              (loop (cdr segments) bestseg best-score)))                       ;; 89-90. bestMatchingSegment = segment; bestScore = nAPS(t-1, segment)
        (values segments bestseg))))                                           ;; 92. return bestMatchingSegment
                                                                                            ;
(define (sorted-segs segs)               ;; (listof Segment) -> (listof Segment)
  ;; produce list of segments sorted by cellx or flatx
  (list-sort (lambda (sega segb)
                (cond [(fx<? (seg-cellx sega) (seg-cellx segb))]
                      [(fx=? (seg-cellx sega) (seg-cellx segb))
                          (fx<? (seg-flatx sega) (seg-flatx segb))]
                      [else #f]))
              segs))
                                                                                            ;
(define (skip-col tm column segments)    ;; TM ColX (listof Segment) -> (listof Segment)
  ;; step to next column's segments
  (if (null? segments) segments
    (let loop ((segments segments))
      (if (fx=? column (segs->colx tm segments))
          (loop (cdr segments))
          segments))))
                                                                                            ;
(define (destroy-segment tm segment)     ;; TM Segment ->
  ;; make segment's flatx available for re-use
  (let* ( (flatx (seg-flatx segment))
          (null-segment (make-seg -1 flatx 0 0 (make-synapses 0))))
    (vector-set! (tm-seg-index tm) flatx null-segment)
    (tm-free-flatx-set! tm (cons flatx (tm-free-flatx tm)))))
                                                                                            ;
(define (cellx->colx tm cellx)           ;; TM CellX -> ColX
  (fxdiv cellx (tm-cells-per-column tm)))
                                                                                            ;
(define (seg-colx tm segment)            ;; TM Segment -> ColX
  (cellx->colx tm (seg-cellx segment)))
                                                                                            ;
(define (colxs->colx tm colxs)           ;; TM (listof ColX) -> ColX
  ;; produce first column index from list, or terminating value if list null
  (if (null? colxs) (tm-num-columns tm) 
      (car colxs)))
                                                                                            ;
(define (segs->cellx tm segments)        ;; TM (listof Segment) -> CellX
  ;; produce cell index for first segment in list, or terminating value if list null
  (if (null? segments) (tm-num-cells tm)
      (seg-cellx (car segments))))
                                                                                            ;
(define (segs->colx tm segments)         ;; TM (listof Segment) -> ColX
  ;; produce column index for first segment in list, or terminating value if list null
  (if (null? segments) (tm-num-columns tm)
      (cellx->colx tm (seg-cellx (car segments)))))
                                                                                            ;
(define (add-to-pre-index                ;; CellX Segment (CellVecOf {FlatX}) ->
          prex seg pre-index)
  ;; add to the list of segment flatxs with synapses from the pre-synaptic cell
  (let ((pre-list (vector-ref pre-index prex)))
    (when (or (null? pre-list) (not (fx=? (seg-flatx seg) (car pre-list))))
      (vector-set! pre-index prex (cons (seg-flatx seg) pre-list)))))
                                                                                            ;
;; --- Synapses ---
                                                                                            ;
(define %syn-least% (least-fixnum))
                                                                                            ;
(define (make-synapse prex perm)         ;; PreX Perm -> Synapse
  (+ (* prex x10k) %syn-least% perm))
                                                                                            ;
(define (syn-prex synapse)               ;; Synapse -> PreX
  (div (- synapse %syn-least%) x10k))
                                                                                            ;
(define (syn-perm synapse)               ;; Synapse -> Perm
  (mod (- synapse %syn-least%) x10k))
                                                                                            ;
  (define synapses:       vector)
  (define make-synapses   make-vector)
  (define synapses-ref    vector-ref)
  (define synapses-set!   vector-set!)
  (define synapses-length vector-length)
  (define list->synapses  list->vector)
  (define synapses->list  vector->list)
                                                                                            ;
(define (build-synapses n proc)          ;; Nat (Nat -> Synapse) -> Synapses
  (let ((synapses (make-synapses n)))
    (do ((i 0 (add1 i))) ((fx=? i n) synapses)
      (synapses-set! synapses i (proc i)))))
                                                                                            ;
(define (in-synapses? cellx synapses)    ;; CellX (vectorof Synapse) -> Boolean
  ;; produce whether cellx equals any pre-synaptic-cell index of synapses
  (let* ( (syn-low  (make-synapse cellx %min-perm%))
          (syn-high (fx+ syn-low %max-perm%)))
    (let search ((left 0) (right (fx- (synapses-length synapses) 1)))
      (if (fx>? left right) #f
        (let ((mid (fxdiv (fx+ left right) 2)))
          (cond 
            [ (fx<? (synapses-ref synapses mid) syn-low ) (search (add1 mid) right) ]
            [ (fx<? syn-high (synapses-ref synapses mid)) (search left (fx- mid 1)) ]
            [ else #t]))))))
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
   
)