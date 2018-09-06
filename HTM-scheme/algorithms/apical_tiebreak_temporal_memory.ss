#!r6rs

;; === HTM-scheme Apical Tiebreak Temporal Memory Copyright 2017 Roger Turner. ===
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Based on code from Numenta Platform for Intelligent Computing (NuPIC) ;;
  ;; which is Copyright (C) 2017, Numenta, Inc.                            ;;
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

  ;; Extends the HTM-scheme Temporal Memory library.
  ;; Translated from numenta htmresearch/.../apical_tiebreak_temporal_memory.py,
  ;; htmresearch/.../numpy_helpers.py, nupic-core/.../SparseMatrixConnections.cpp --
  ;; see comments there for descriptions of functions and parameters.
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.

(library (HTM-scheme HTM-scheme algorithms apical_tiebreak_temporal_memory)
                                                                                            ;
(export
  ;; for tests in sequence_memory
  create-segment
  grow-synapses
  tm-basal-connections
  tm-basal-pre-index
  syn-perm
  syn-prex
  seg-synapses
  seg-synapses-set!
  tm-next-flatx
  tm-predicted-cells
  
  map-segments-to-cells
  tm-active-basal-segments
  tm-active-apical-segments

  tm:permanence
  tm:synapse
  tm:segment
  tm:get-active-cols
  (rename
    (tm                   tm:tm)
    (make-tm              tm:construct)
    (reset                tm:reset)
    (activate-cells       tm:activate-cells)
    (depolarize-cells     tm:depolarize-cells)
    (get-active-cells     tm:get-active-cells)
    (get-winner-cells     tm:get-winner-cells)))
    
                                                                                              ;
(import
  (rnrs)
  (HTM-scheme HTM-scheme algorithms htm_prelude)
  (rename (HTM-scheme HTM-scheme algorithms temporal_memory)
    (tm-active-segments             tm-active-basal-segments)
    (tm-active-segments-set!        tm-active-basal-segments-set!)
    (tm-matching-segments           tm-matching-basal-segments)
    (tm-matching-segments-set!      tm-matching-basal-segments-set!)
    (tm-cells                       tm-basal-connections)
    (tm-pre-index                   tm-basal-pre-index)
    (tm-predicted-segment-decrement tm-basal-predicted-segment-decrement)))
    
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
;; CellX        = Nat [0 - (div (expt 2 (fixnum-width)) 10000)] cell index
;; CellVecOf    = Vector indexed by CellX
;; ColX         = Nat [0 - MAXCOL], column index of cell
;; Connections  = (CellVecOf (listof Segment)), basal and apical segments for cells
;; TM           = Record: tm parameters, cells

;; === Parameters, Data, Convenience Functions ===
                                                                                            ;
(define-record-type tm                     ;; TM
  (parent tm:tm)
  (fields
    basal-input-size                       ;; Nat
    apical-input-size                      ;; Nat
    ;;cells-per-column                     ;; Nat [inherited from tm:tm]
    reduced-basal-threshold                ;; Nat
    sample-size                            ;; Nat
    ;;basal-predicted-segment-decrement    ;; Permanence [renamed tm:tm predicted-segment-decrement]
    apical-predicted-segment-decrement     ;; Permanence
    (mutable predicted-cells)              ;; (listof CellX)
    (mutable predicted-active-cells)       ;; (listof CellX) 
    ;;(mutable active-basal-segments)      ;; (listof Segment) [renamed tm:tm active-segments]
    ;;(mutable matching-basal-segments)    ;; (listof Segment) [renamed tm:tm matching-segments]
    ;;basal-connections                    ;; (CellVecOf (listof Segment)) [renamed tm cells]
    (mutable basal-potential-overlaps)     ;; (vectorof Nat) [indexed by flatx]
    (mutable active-apical-segments)       ;; (listof Segment)
    (mutable matching-apical-segments)     ;; (listof Segment)
    apical-connections                     ;; Connections
    apical-pre-index                       ;; (CellVecOf {FlatX})
    (mutable apical-potential-overlaps)    ;; (vectorof Nat) [indexed by flatx]
    (mutable use-apical-tiebreak)          ;; Boolean
    (mutable use-apical-modulation-basal-threshold)) ;; Boolean
    
  (protocol
    (lambda (pargs->new)                   ;; Nat Nat (listof KWarg) -> TM
      (lambda (column-count cells-per-column . kwargs)
        (let ((column-dimensions (list column-count))
              (num-cells         (* column-count cells-per-column)))
          (apply (apply pargs->new (append (list column-dimensions cells-per-column) kwargs))
                 (key-word-args
                    (append `(
                              [apical-connections     . ,(make-vector num-cells '())]
                              [apical-pre-index       . ,(make-vector num-cells '())]
                              )
                      kwargs)
                    attm-defaults #f)))))))
                                                                                            ;
(define attm-defaults                      ;; (listof KWarg)
  `(
    [basal-input-size                   . 0]
    [apical-input-size                  . 0]
    [reduced-basal-threshold            . 13]
    [sample-size                        . 20]
    [apical-predicted-segment-decrement . ,(tm:permanence 0.0)]
    [predicted-cells                    . ()]
    [predicted-active-cells             . ()]
    [basal-potential-overlaps           . #()]
    [active-apical-segments             . ()]
    [matching-apical-segments           . ()]
    [apical-connections                 . #()]
    [apical-pre-index                   . #()]
    [apical-potential-overlaps          . #()]
    [use-apical-tiebreak                . #t]
    [use-apical-modulation-basal-threshold . #t]))
                                                                                            ;
(define (in1d x1s x2s)                     ;; {X} {X} -> Bitwise
  ;; produce index mask of x1s that are also in x2s
  (fold-left 
    (lambda (acc e1 e1x)
      (if (memv e1 x2s)
        (bitwise-copy-bit acc e1x 1)
        acc))
    0 x1s (build-list (length x1s) id)))
                                                                                            ;
(define (include-by-mask xs mask)          ;; {X} Bitwise -> {X}
  ;; extract elements of xs corresponding to 1 bits in mask
  (vector->list
    (vector-refs 
      (list->vector xs) 
      (list->vector (bitwise->list mask)))))
                                                                                            ;
(define (exclude-by-mask xs mask)          ;; {X} Bitwise -> {X}
  (vector->list
    (vector-refs 
      (list->vector xs) 
      (list->vector 
        (bitwise->list
          (bitwise-bit-field (bitwise-not mask) 0 (length xs)))))))
                                                                                            ;
(define (intersect1d l1 l2)                ;; {Nat} {Nat} -> {Nat}
  ;; produce intersection of sorted lists of naturals, ie
  ;; (bitwise->list (bitwise-and (list->bitwise l1) (list->bitwise l2)))
  (let loop ((l1 l1) (l2 l2) (result (list)))
    (cond [(or (null? l1) (null? l2)) result]
          [(fx=? (car l1) (car l2))
              (loop (cdr l1) (cdr l2) (cons (car l1) result))]
          [(fx<? (car l1) (car l2))
              (loop (cdr l1) l2 result)]
          [else (loop l1 (cdr l2) result)])))
                                                                                            ;
(define (setdiff1d l1 l2)                  ;; {Nat} {Nat} -> {Nat}
  ;; produce difference of sorted lists of naturals, ie
  ;; (bitwise->list (bitwise-and (list->bitwise l1) (bitwise-not (list->bitwise l2))))
  (let loop ((l1 l1) (l2 l2) (result (list)))
    (cond [(null? l1) result]
          [(null? l2) (append result l1)]
          [(fx=? (car l1) (car l2))
              (loop (cdr l1) (cdr l2) result)]
          [(fx<? (car l1) (car l2))
              (loop (cdr l1) l2 (cons (car l1) result))]
          [else (loop l1 (cdr l2) (cons (car l1) result))])))
                                                                                            ;
(define (map-segments-to-cells segments)   ;; {Segment} -> {CellX}
  (map seg-cellx segments))
                                                                                            ;
(define (filter-segments-by-cell           ;; {Seg} {CellX} -> {Seg}
          segments cellxs)
  ;; Return the subset of segments that are on the provided cells.
  (filter
    (lambda (segment)
      (memv (seg-cellx segment) cellxs))
    segments))
                                                                                            ;
(define (map-segments-to-synapse-counts    ;; {Seg} -> {Nat}
          segments)
  (map
    (lambda (segment)
      (synapses-length (seg-synapses segment)))
    segments))
                                                                                            ;
(define (/cpc tm)                          ;; TM -> (CellX -> ColX)
  (lambda (cellx) 
    (fxdiv cellx (tm-cells-per-column tm))))
                                                                                            ;
(define (set-compare tm a b)               ;; TM {CellX} {ColX} -> {CellX} {CellX}
  ;; produce intersection (a&b), difference (b-a), using (a/cpc) for compare
  (let ((a-within-b-mask (in1d (map (/cpc tm) a) b))
        (b-within-a-mask (in1d b (map (/cpc tm) a))))
    (values (include-by-mask a a-within-b-mask) (exclude-by-mask b b-within-a-mask))))
                                                                                           ;
(define (argmax-multi a group-keys)        ;; {Number} {Number} -> {Nat}
  ;; gets the indices of the max values of each group in 'a', grouping the
  ;; elements by their corresponding value in 'groupKeys'.
  (let per-group ((a a) (group-keys group-keys) (index 0) (result '()))
    (cond
      [ (null? group-keys) (reverse result) ]
      [ else
        (let ((group (car group-keys)))
          (let next ((a a) (group-keys group-keys) (index index) (max-in-group (least-fixnum)) (index-of-max #f))
            (cond
              [ (null? group-keys) (per-group a group-keys index (cons index-of-max result)) ]
              [ (fx=? (car group-keys) group)
                  (if (fx>? (car a) max-in-group)
                    (next (cdr a) (cdr group-keys) (add1 index) (car a) index)
                    (next (cdr a) (cdr group-keys) (add1 index) max-in-group index-of-max)) ]
              [ else (per-group a group-keys index (cons index-of-max result)) ]))) ])))
                                                                                            ;
(define (get-all-cells-in-columns tm colxs);; TM {ColX} -> {CellX}
  ;; Calculate all cell indices in the specified columns.
  (let ((cpc (tm-cells-per-column tm)))
    (apply append 
      (map
        (lambda (colx)
          (let ((first (fx* colx cpc)))
            (build-list cpc (lambda (i) (fx+ first i)))))
        colxs))))
                                                                                            ;
(define (adjust-synapses tm                ;; TM Connections {Seg} (vectorof CellX) ->
          connections segments active-input)
  ;; For each specified segment, update the permanence of each synapse
  ;; according to whether the synapse would be active given the specified
  ;; active inputs.
  (for-each
    (lambda (segment)
      (adapt-segment tm segment (list->vector active-input) connections))
    segments))
                                                                                            ;
(define (adjust-active-synapses tm         ;; TM Connections {Seg} {CellX} Permanence ->
          connections segments active-input permanence-delta)
  ;; For each specified segment, add a delta to the permanences of the
  ;; synapses that would be active given the specified active inputs.
  (for-each
    (lambda (segment)
      (adapt-segment tm segment (list->vector active-input) connections
                     (lambda (p)
                       (clip-max (clip-min (fx+ p permanence-delta))))
                     0))
    segments))
                                                                                            ;
(define (grow-synapses-to-sample tm        ;; TM {Seg} {CellX} Nat|{Nat} (CellVecOf {FlatX}) ->
          segments inputs sample-size pre-index)
  ;; For each segment, grow synapses to a random subset of the
  ;; inputs that aren't already connected to the segment.
  (if (fixnum? sample-size)
    (for-each
      (lambda (segment)
        (grow-synapses tm segment sample-size inputs pre-index))
      segments)
    (for-each
      (lambda (segment sample-size)
        (grow-synapses tm segment sample-size inputs pre-index))
      segments sample-size)))
                                                                                            ;
(define (grow-synapses                   ;; TM Segment Nat (listof CellX) [(CellVecOf {FlatX})] ->
          tm segment n-desired-new-synapses inputs pre-index)
  ;; grow synapses to all inputs that aren't already connected to the segment
  (let* ( (synapses     (seg-synapses segment))
          (num-synapses (synapses-length synapses))
          (candidates   (let loop ((cs '()) (is inputs))
                          (cond [ (null? is) cs ]
                                [ (in-synapses? (car is) synapses)
                                    (loop cs (cdr is)) ]
                                [ else (loop (cons (car is) cs) (cdr is)) ])))
          (n-actual     (fxmin n-desired-new-synapses (length candidates))))
      (when (positive? n-actual)
        (let ((new-prexs (vector-sample (list->vector candidates) n-actual)))
          (seg-synapses-set! segment
            (build-synapses (fx+ num-synapses n-actual)
              (lambda (sx)
                ;; copy retained synapses, make n-actual new ones
                (if (fx<? sx num-synapses)
                    (synapses-ref synapses sx)
                    (let ((new-prex (vector-ref new-prexs (fx- sx num-synapses))))
                      (add-to-pre-index new-prex segment pre-index)
                      (tm:synapse new-prex (tm-initial-permanence tm)))))))
          ;; keep synapses sorted for binary search in compute-activity
          (vector-sort! fx<? (seg-synapses segment))))))
                                                                                            ;
(define (create-segments tm                ;; TM Connections {CellX} -> {Seg}
          connections cellxs)
  (map
    (lambda (cellx)
      (create-segment tm cellx connections))
    cellxs))

;; === Apical Tiebreak Temporal Memory Algorithm ===
                                                                                            ;
(define (reset tm)                         ;; TM ->
  ;; Clear all cell and segment activity.
  (tm:reset tm)
  (tm-predicted-cells-set! tm                '())
  (tm-predicted-active-cells-set! tm         '())
  (tm-active-apical-segments-set! tm         '())
  (tm-matching-apical-segments-set! tm       '())
  (tm-basal-potential-overlaps-set! tm  (vector))
  (tm-apical-potential-overlaps-set! tm (vector)))
                                                                                            ;
(define (depolarize-cells                  ;; TM {CellX} {CellX} Boolean ->
          tm basal-input apical-input learn)
  ;; Calculate predictions.
  (let-values ([(active-apical-segments matching-apical-segments apical-potential-overlaps)
                (calculate-apical-segment-activity 
                    tm (list->vector apical-input) (tm-apical-pre-index tm))])
    (let ((reduced-basal-threshold-cells
            (if (or learn (tm-use-apical-modulation-basal-threshold tm))
              '()
              (map-segments-to-cells active-apical-segments))))
      (let-values ([(active-basal-segments matching-basal-segments basal-potential-overlaps)
                    (calculate-basal-segment-activity tm (list->vector basal-input) 
                        (tm-basal-pre-index tm) reduced-basal-threshold-cells)])
        (tm-predicted-cells-set! tm 
          (calculate-predicted-cells tm active-basal-segments active-apical-segments))
        (tm-active-basal-segments-set!     tm active-basal-segments)
        (tm-active-apical-segments-set!    tm active-apical-segments)
        (tm-matching-basal-segments-set!   tm (sorted-segs matching-basal-segments))
        (tm-matching-apical-segments-set!  tm (sorted-segs matching-apical-segments))
        (tm-basal-potential-overlaps-set!  tm basal-potential-overlaps)
        (tm-apical-potential-overlaps-set! tm apical-potential-overlaps)))))
                                                                                            ;
(define (activate-cells                    ;; TM {ColX} {CellX} {CellX} {CellX} {CellX} Boolean ->
          tm active-columns basal-reinforce-candidates apical-reinforce-candidates
          basal-growth-candidates apical-growth-candidates learn)
  ;; Activate cells in the specified columns, using the result of the previous
  ;; 'depolarizeCells' as predictions. Then learn.
  (let*-values (
    ;; Calculate active cells
    [ ( correct-predicted-cells bursting-columns)
          (set-compare tm (tm-predicted-cells tm) active-columns)]
    [ ( new-active-cells) 
          (append correct-predicted-cells (get-all-cells-in-columns tm bursting-columns)) ]
    ;; Calculate learning
    [ ( learning-active-basal-segments
        learning-matching-basal-segments
        basal-segments-to-punish
        new-basal-segment-cells
        learning-cells)
          (calculate-basal-learning tm
            active-columns bursting-columns correct-predicted-cells
            (tm-active-basal-segments tm) (tm-matching-basal-segments tm)
            (tm-basal-potential-overlaps tm)) ]
    [ ( learning-active-apical-segments
        learning-matching-apical-segments
        apical-segments-to-punish
        new-apical-segment-cells)
          (calculate-apical-learning tm
            learning-cells active-columns
            (tm-active-apical-segments tm) (tm-matching-apical-segments tm)
            (tm-apical-potential-overlaps tm)) ] )
    (when learn
      ;; Learn on existing segments
      (learn-on-segments tm (tm-basal-connections tm) learning-active-basal-segments
                         basal-reinforce-candidates basal-growth-candidates
                         (tm-basal-potential-overlaps tm) (tm-basal-pre-index tm))
      (learn-on-segments tm (tm-basal-connections tm) learning-matching-basal-segments
                         basal-reinforce-candidates basal-growth-candidates
                         (tm-basal-potential-overlaps tm) (tm-basal-pre-index tm))
      (learn-on-segments tm (tm-apical-connections tm) learning-active-apical-segments
                         apical-reinforce-candidates apical-growth-candidates
                         (tm-apical-potential-overlaps tm) (tm-apical-pre-index tm))
      (learn-on-segments tm (tm-apical-connections tm) learning-matching-apical-segments
                         apical-reinforce-candidates apical-growth-candidates
                         (tm-apical-potential-overlaps tm) (tm-apical-pre-index tm))
      ;; Punish incorrect predictions
      (unless (zero? (tm-basal-predicted-segment-decrement tm))
        (adjust-active-synapses tm (tm-basal-connections tm)
                                basal-segments-to-punish basal-reinforce-candidates
                                (- (tm-basal-predicted-segment-decrement tm))))
      (unless (zero? (tm-apical-predicted-segment-decrement tm))
        (adjust-active-synapses tm (tm-apical-connections tm)
                                apical-segments-to-punish apical-reinforce-candidates
                                (- (tm-apical-predicted-segment-decrement tm))))
      ;; Grow new segments
      (when (positive? (length basal-growth-candidates))
        (learn-on-new-segments tm (tm-basal-connections tm) new-basal-segment-cells 
                               basal-growth-candidates (tm-basal-pre-index tm)))
      (when (positive? (length apical-growth-candidates))
        (learn-on-new-segments tm (tm-apical-connections tm) new-apical-segment-cells
                               apical-growth-candidates (tm-apical-pre-index tm))))
      ;; Save the results
    (tm-active-cells-set! tm (list->vector new-active-cells))
    (vector-sort! fx<? (tm-active-cells tm))
    (tm-winner-cells-set! tm (list-sort fx<? learning-cells))
    (tm-predicted-active-cells-set! tm correct-predicted-cells)))
                                                                                          ;
(define (calculate-basal-learning tm       ;; TM {ColX} {ColX} {CellX} {Seg} {Seg} {Nat}
          active-columns bursting-columns correct-predicted-cells
          active-basal-segments matching-basal-segments basal-potential-overlaps)
  ;; Basic Temporal Memory learning. Correctly predicted cells always have
  ;; active basal segments, and we learn on these segments. In bursting
  ;; columns, we either learn on an existing basal segment, or we grow a new one.
  ;; The only influence apical dendrites have on basal learning is: the apical
  ;; dendrites influence which cells are considered "predicted". So an active
  ;; apical dendrite can prevent some basal segments in active columns from
  ;; learning.
  ;; Correctly predicted columns
  (let* ( (learning-active-basal-segments 
            (filter-segments-by-cell active-basal-segments correct-predicted-cells))
          (cells-for-matching-basal
            (map-segments-to-cells matching-basal-segments))
          (matching-cells (unique fx=? cells-for-matching-basal)))
    (let-values ([(matching-cells-in-bursting-columns bursting-columns-with-no-match)
                  (set-compare tm matching-cells bursting-columns)])
      (let* ( (learning-matching-basal-segments
                (choose-best-segment-per-column tm matching-cells-in-bursting-columns
                  matching-basal-segments basal-potential-overlaps))
              (new-basal-segment-cells
                (get-cells-with-fewest-segments tm bursting-columns-with-no-match))
              (learning-cells
                (append correct-predicted-cells
                        (map-segments-to-cells learning-matching-basal-segments)
                        new-basal-segment-cells))
              ;; Incorrectly predicted columns
              (correct-matching-basal-mask
                (in1d (map (/cpc tm) cells-for-matching-basal) active-columns))
              (basal-segments-to-punish
                (exclude-by-mask matching-basal-segments correct-matching-basal-mask)))
        (values
          learning-active-basal-segments
          learning-matching-basal-segments
          basal-segments-to-punish
          new-basal-segment-cells
          learning-cells)))))
                                                                                          ;
(define (calculate-apical-learning tm      ;; TM {CellX} {ColX} {Seg} {Seg} {Nat}
          learning-cells active-columns  
          active-apical-segments matching-apical-segments apical-potential-overlaps)
  ;; Calculate apical learning for each learning cell.
  ;; The set of learning cells was determined completely from basal segments.
  ;; Do all apical learning on the same cells.
  ;; Learn on any active segments on learning cells. For cells without active
  ;; segments, learn on the best matching segment. For cells without a matching
  ;; segment, grow a new segment.
  (let* ( (learning-active-apical-segments 
            (filter-segments-by-cell active-apical-segments learning-cells))
          (learning-cells-without-active-apical
            (setdiff1d learning-cells (map-segments-to-cells learning-active-apical-segments)))
          (cells-for-matching-apical
            (map-segments-to-cells matching-apical-segments))
          (learning-cells-with-matching-apical
            (intersect1d learning-cells-without-active-apical cells-for-matching-apical))
          (learning-matching-apical-segments
            (choose-best-segment-per-cell learning-cells-with-matching-apical
                                          matching-apical-segments apical-potential-overlaps))
          (new-apical-segment-cells
            (setdiff1d learning-cells-without-active-apical learning-cells-with-matching-apical))
          (correct-matching-apical-mask
            (in1d (map (/cpc tm) cells-for-matching-apical) active-columns))
          (apical-segments-to-punish
            (exclude-by-mask matching-apical-segments correct-matching-apical-mask)))
        (values
          learning-active-apical-segments
          learning-matching-apical-segments
          apical-segments-to-punish
          new-apical-segment-cells)))
                                                                                          ;
(define (calculate-apical-segment-activity ;; TM (vectorof CellX) (CellVecOf {FlatX}) -> {Segment} {Segment} (vectorof Nat)
          tm active-input pre-index)
  ;; Calculate the active and matching apical segments for this timestep.
  (compute-activity tm active-input pre-index 0 '()))
                                                                                            ;
(define (calculate-basal-segment-activity  ;; TM (vectorof CellX) (CellVecOf {FlatX}) {CellX} -> {Segment} {Segment} (vectorof Nat)
          tm active-input pre-index reduced-basal-threshold-cells)
  ;; Calculate the active and matching basal segments for this timestep.
  (compute-activity tm active-input pre-index 
                    (tm-reduced-basal-threshold tm) reduced-basal-threshold-cells))
                                                                                            ;
(define (calculate-predicted-cells         ;; TM {Segment} {Segment} -> {CellX}
          tm active-basal-segments active-apical-segments)
  ;; Calculate the predicted cells, given the set of active segments.
  (let* ( (cells-for-basal-segments  (map-segments-to-cells active-basal-segments))
          (cells-for-apical-segments (map-segments-to-cells active-apical-segments))
          (fully-depolarized-cells         ;; cells with both types of segments active
            (intersect1d cells-for-basal-segments cells-for-apical-segments))
          (partly-depolarized-cells        ;; cells with basal only
            (setdiff1d cells-for-basal-segments fully-depolarized-cells))
          (inhibited-mask
            (in1d 
              (map (/cpc tm) partly-depolarized-cells)
              (map (/cpc tm) fully-depolarized-cells)))
          (predicted-cells
            (list-sort fx<?
              (append fully-depolarized-cells
                      (exclude-by-mask partly-depolarized-cells inhibited-mask)))))
    predicted-cells))
                                                                                            ;
(define (learn-on-segments tm              ;; TM Connections {Seg} (vectorof CellX) {CellX} (vectorof Nat) (CellVecOf {FlatX}) ->
          connections learning-segments active-input growth-candidates 
          potential-overlaps pre-index)
  ;; Adjust synapse permanences, grow new synapses, and grow new segments.
  (adjust-synapses tm connections learning-segments active-input)
    ;; Grow new synapses. Calculate "maxNew", the maximum number of synapses to
    ;; grow per segment. "maxNew" might be a number or it might be a list of numbers.
  (let ((max-new
          (if (fx=? -1 (tm-sample-size tm))
            (length growth-candidates)
            (map
              (lambda (segment)
                (fx- (tm-sample-size tm) (vector-ref potential-overlaps (seg-flatx segment))))
              learning-segments))))
    (let ((max-new
            (if (fx=? -1 (tm-max-synapses-per-segment tm))
              max-new
              (let* ((synapse-counts
                       (map-segments-to-synapse-counts learning-segments))
                     (num-synapses-to-reach-max
                       (map fx- 
                          (make-list (length synapse-counts) (tm-max-synapses-per-segment tm))
                          synapse-counts)))
                (map
                  (lambda (mn nstrm)
                    (if (fx<=? mn nstrm)
                      mn
                      nstrm))
                  max-new num-synapses-to-reach-max)))))
      (grow-synapses-to-sample tm learning-segments growth-candidates max-new pre-index))))
                                                                                            ;
(define (learn-on-new-segments tm          ;; TM Connections {CellX} {CellX} (CellVecOf {FlatX}) ->
          connections new-segment-cells growth-candidates pre-index)
  (let* ( (num-new-synapses (length growth-candidates))
          (num-new-synapses (if (fx=? -1 (tm-sample-size tm))
                              num-new-synapses
                              (fxmin num-new-synapses (tm-sample-size tm))))
          (num-new-synapses (if (fx=? -1 (tm-max-synapses-per-segment tm))
                              num-new-synapses
                              (fxmin num-new-synapses (tm-max-synapses-per-segment tm))))
          (new-segments (create-segments tm connections new-segment-cells)))
    (grow-synapses-to-sample tm new-segments growth-candidates num-new-synapses pre-index)))
                                                                                            ;
(define (choose-best-segment-per-cell      ;; {CellX} {Seg} (vectorof Nat) -> {Seg}
          cells all-matching-segments potential-overlaps)
  ;; For each specified cell, choose its matching segment with largest number
  ;; of active potential synapses. When there's a tie, the first segment wins.
  (let next-cell ((cells cells) (learning-segments '()))
    (if (null? cells) 
      (reverse learning-segments)
      (let next-segment ( (segments all-matching-segments) 
                          (best-seg #f) 
                          (max-overlap (least-fixnum)))
        (cond
          [ (null? segments) (next-cell (cdr cells) (cons best-seg learning-segments)) ]
          [ (fx=? (car cells) (seg-cellx (car segments)))
              (let ((overlap (vector-ref potential-overlaps (seg-flatx (car segments)))))
                (if (fx>? overlap max-overlap)
                    (next-segment (cdr segments) (car segments) overlap)
                    (next-segment (cdr segments) best-seg max-overlap))) ]
          [ else (next-segment (cdr segments) best-seg max-overlap) ] )))))
                                                                                            ;
(define (choose-best-segment-per-column tm ;; TM {CellX} {Seg} (vectorof Nat) -> {Seg}
          matching-cells all-matching-segments potential-overlaps)
  ;; For all the columns covered by 'matchingCells', choose the column's matching
  ;; segment with largest number of active potential synapses. When there's a
  ;; tie, the first segment wins.
  (let* ( (candidate-segments
            (filter-segments-by-cell all-matching-segments matching-cells))
          (cell-scores
            (map
              (lambda (segment)
                (vector-ref potential-overlaps (seg-flatx segment)))
              candidate-segments))
          (columns-for-candidates
            (map (/cpc tm) (map-segments-to-cells candidate-segments)))
          (one-per-column-filter 
            (argmax-multi cell-scores columns-for-candidates))
          (learning-segments
            (vector->list (vector-refs (list->vector candidate-segments)
                                       (list->vector one-per-column-filter)))))
    learning-segments))
                                                                                            ;
(define (get-cells-with-fewest-segments tm ;; TM {ColX} -> {CellX}
          columnxs)
  ;; For each column, get the cell that has the fewest total basal segments.
  ;; Break ties randomly.
  (map
    (lambda (colx)
      (let loop ((cellx (fx* colx (tm-cells-per-column tm)))
                 (candidate-cells '())
                 (fewest-segs (greatest-fixnum)))
        (if (fx<? cellx (fx* (add1 colx) (tm-cells-per-column tm)))
          (let ((n-segs (length (vector-ref (tm-basal-connections tm) cellx))))
            (cond
              [ (fx<? n-segs fewest-segs)
                  (loop (add1 cellx) (list cellx) n-segs) ]
              [ (fx=? n-segs fewest-segs)
                  (loop (add1 cellx) (cons cellx candidate-cells) n-segs) ]
              [ else (loop (add1 cellx) candidate-cells fewest-segs) ] ))
          (list-ref candidate-cells (random (length candidate-cells))))))
    columnxs))
                                                                                            ;
(define (get-active-cells tm)              ;; TM -> {CellX}
  (vector->list (tm-active-cells tm)))
                                                                                            ;
(define (get-predicted-active-cells tm)    ;; TM -> {CellX}
  (tm-predicted-active-cells tm))
                                                                                            ;
(define (get-winner-cells tm)              ;; TM -> {CellX}
  (tm-winner-cells tm))
                                                                                            ;
)