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

  ;; Translated from numenta htmresearch/.../apical_tiebreak_temporal_memory.py,
  ;; htmresearch/.../numpy_helpers.py, nupic-core/.../SparseMatrixConnections.cpp --
  ;; see comments there for descriptions of functions and parameters.
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.

(library (HTM-scheme HTM-scheme algorithms apical_tiebreak_temporal_memory)
;(library (apical_tiebreak_temporal_memory)
                                                                                            ;
(export
  tm
  make-tm
  reset
  activate-cells
  depolarize-cells
  get-active-cells
  get-predicted-active-cells
  get-winner-cells
  get-predicted-cells
  map-segments-to-cells
  get-active-basal-segments
  get-active-apical-segments
  (rename                                ;; for temporal_memory_test
    (create-segment       attm:create-segment)
    (tm-basal-connections attm:tm-basal-connections)
    (grow-synapses        attm:grow-synapses)
    (tm-basal-pre-index   attm:tm-basal-pre-index)
    (tm-next-flatx        attm:tm-next-flatx)
    ))
                                                                                            ;
(import
  (rnrs)
  (except (HTM-scheme HTM-scheme algorithms htm_prelude) add1 search)
  (HTM-scheme HTM-scheme algorithms htm_concept))
  
(define (add1 n)                         ;; Fixnum -> Fixnum
  (fx+ 1 n))
                                                                                            ;
(define (search synapses syn-low syn-high)
  (let search ((left 0) (right (fx- (vector-length synapses) 1)))
    (if (fx>? left right) #f
      (let* ( (mid (fxdiv (fx+ left right) 2))
              (synapse (vector-ref synapses mid))) 
        (cond 
          [ (fx<? synapse  syn-low) (search (add1 mid) right) ]
          [ (fx<? syn-high synapse) (search left (fx- mid 1)) ]
          [else synapse])))))

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
(define-record-type tm                   ;; TM
  (fields
    column-count                         ;; Nat  The number of minicolumns
    basal-input-size                     ;; Nat  The number of bits in the basal input
    apical-input-size                    ;; Nat  The number of bits in the apical input
    cells-per-column                     ;; Nat 
    activation-threshold                 ;; Nat
    reduced-basal-threshold              ;; Nat
    initial-permanence                   ;; Perm
    connected-permanence                 ;; Perm
    min-threshold                        ;; Nat
    sample-size                          ;; Nat
    permanence-increment                 ;; Perm
    permanence-decrement                 ;; Perm
    basal-predicted-segment-decrement    ;; Perm
    apical-predicted-segment-decrement   ;; Perm
    max-synapses-per-segment             ;; Fixnum [-1 => unlimited]
    seed                                 ;; Integer
    use-apical-tiebreak                  ;; Boolean
    use-apical-modulation-basal-threshold ;; Boolean
    (mutable basal-connections)          ;; (CellVecOf (listof Segment))
    (mutable apical-connections)         ;; (CellVecOf (listof Segment))
    (mutable active-cells)               ;; (vectorof CellX)
    (mutable winner-cells)               ;; (listof CellX)
    (mutable predicted-cells)            ;; (listof CellX)
    (mutable predicted-active-cells)     ;; (listof CellX) 
    (mutable active-basal-segments)      ;; (listof Segment)
    (mutable active-apical-segments)     ;; (listof Segment)
    (mutable matching-basal-segments)    ;; (listof Segment)
    (mutable matching-apical-segments)   ;; (listof Segment)
    (mutable basal-potential-overlaps)   ;; (vectorof Nat) [indexed by flatx]
    (mutable apical-potential-overlaps)  ;; (vectorof Nat) [indexed by flatx]
    (mutable basal-pre-index)            ;; (CellVecOf (listof FlatX))
    (mutable apical-pre-index)           ;; (CellVecOf (listof FlatX))
    (mutable seg-index)                  ;; (vectorof Segment) [indexed by flatx]
    (mutable next-flatx)                 ;; Nat: next available index in seg-index
    (mutable free-flatx))                ;; (listof Nat): indices available for re-use
  (protocol
    (lambda (new)                        ;; (listof KWarg) -> TM
      (lambda (kwargs)
        (let* ( (tm (apply new (key-word-args kwargs attm-defaults)))
                (num-cells (* (tm-column-count tm) (tm-cells-per-column tm))))
          (tm-basal-connections-set!  tm (make-vector num-cells '()))
          (tm-apical-connections-set! tm (make-vector num-cells '()))
          (tm-basal-pre-index-set!    tm (make-vector num-cells '()))
          (tm-apical-pre-index-set!   tm (make-vector num-cells '()))
          (tm-seg-index-set!          tm (make-vector (tm-column-count tm)))
          tm)))))
                                                                                            ;
(define attm-defaults `(                 ;; (listof KWarg)
  [column-count                          . 2048]
  [basal-input-size                      . 0]
  [apical-input-size                     . 0]
  [cells-per-column                      . 32]
  [activation-threshold                  . 13]
  [reduced-basal-threshold               . 13]
  [initial-permanence                    . ,(perm 0.21)]
  [connected-permanence                  . ,(perm 0.50)]
  [min-threshold                         . 10]
  [sample-size                           . 20]
  [permanence-increment                  . ,(perm 0.1)]
  [permanence-decrement                  . ,(perm 0.1)]
  [basal-predicted-segment-decrement     . ,(perm 0.0)]
  [apical-predicted-segment-decrement    . ,(perm 0.0)]
  [max-synapses-per-segment              . -1]
  [seed                                  . 42]
  [use-apical-tiebreak                   . #t]
  [use-apical-modulation-basal-threshold . #t]
  [basal-connections                     . #()]
  [apical-connections                    . #()]
  [active-cells                          . #()]
  [winner-cells                          . ()]
  [predicted-cells                       . ()]
  [predicted-active-cells                . ()]
  [active-basal-segments                 . ()]
  [active-apical-segments                . ()]
  [matching-basal-segments               . ()]
  [matching-apical-segments              . ()]
  [basal-potential-overlaps              . #()]
  [apical-potential-overlaps             . #()]
  [basal-pre-index                       . #()]
  [apical-pre-index                      . #()]
  [seg-index                             . #()]
  [next-flatx                            . 0]
  [free-flatx                            . ()]))

;; === Apical Tiebreak Temporal Memory Algorithm ===
                                                                                            ;
(define (reset tm)                       ;; TM ->
  ;; Clear all cell and segment activity.
  (tm-active-cells-set!              tm (vector))
  (tm-winner-cells-set!              tm      '())
  (tm-predicted-cells-set!           tm      '())
  (tm-predicted-active-cells-set!    tm      '())
  (tm-active-basal-segments-set!     tm      '())
  (tm-active-apical-segments-set!    tm      '())
  (tm-matching-basal-segments-set!   tm      '())
  (tm-matching-apical-segments-set!  tm      '())
  (tm-basal-potential-overlaps-set!  tm (vector))
  (tm-apical-potential-overlaps-set! tm (vector)))
                                                                                            ;
(define (depolarize-cells                ;; TM {CellX} {CellX} Boolean ->
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
(define (activate-cells                  ;; TM {ColX} {CellX} {CellX} {CellX} {CellX} Boolean ->
          tm active-columns basal-reinforce-candidates apical-reinforce-candidates
          basal-growth-candidates apical-growth-candidates learning)
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
    (when learning
      ;; Learn on existing segments
      (learn tm (tm-basal-connections tm) learning-active-basal-segments
        basal-reinforce-candidates basal-growth-candidates
        (tm-basal-potential-overlaps tm) (tm-basal-pre-index tm))
      (learn tm (tm-basal-connections tm) learning-matching-basal-segments
        basal-reinforce-candidates basal-growth-candidates
        (tm-basal-potential-overlaps tm) (tm-basal-pre-index tm))
      (learn tm (tm-apical-connections tm) learning-active-apical-segments
        apical-reinforce-candidates apical-growth-candidates
        (tm-apical-potential-overlaps tm) (tm-apical-pre-index tm))
      (learn tm (tm-apical-connections tm) learning-matching-apical-segments
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
    (tm-winner-cells-set! tm learning-cells)
    (tm-predicted-active-cells-set! tm correct-predicted-cells)))
                                                                                          ;
(define (calculate-basal-learning tm     ;; TM {ColX} {ColX} {CellX} {Seg} {Seg} {Nat}
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
              (learning-cells (list-sort fx<?
                (append correct-predicted-cells
                        (map-segments-to-cells learning-matching-basal-segments)
                        new-basal-segment-cells)))
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
(define (calculate-apical-learning tm    ;; TM {CellX} {ColX} {Seg} {Seg} {Nat}
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
(define (calculate-predicted-cells tm    ;; TM {Segment} {Segment} -> {CellX}
          active-basal-segments active-apical-segments)
  ;; Calculate the predicted cells, given the set of active segments.
  (let* ( (cells-for-basal-segments  (map-segments-to-cells active-basal-segments))
          (cells-for-apical-segments (map-segments-to-cells active-apical-segments))
          (fully-depolarized-cells       ;; cells with both types of segments active
            (intersect1d cells-for-basal-segments cells-for-apical-segments))
          (partly-depolarized-cells      ;; cells with basal only
            (setdiff1d cells-for-basal-segments fully-depolarized-cells))
          (inhibited-mask
            (in1d 
              (map (/cpc tm) partly-depolarized-cells)
              (map (/cpc tm) fully-depolarized-cells)))
          (predicted-cells
            (if (tm-use-apical-tiebreak tm)
              (list-sort fx<?
                (append fully-depolarized-cells
                        (exclude-by-mask partly-depolarized-cells inhibited-mask)))
                cells-for-basal-segments)))
    predicted-cells))
                                                                                            ;
(define (learn tm connections            ;; TM Connections {Seg} (vectorof CellX) {CellX} (vectorof Nat) (CellVecOf {FlatX}) ->
          learning-segments active-input growth-candidates potential-overlaps pre-index)
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
(define (learn-on-new-segments tm        ;; TM Connections {CellX} {CellX} (CellVecOf {FlatX}) ->
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
(define (choose-best-segment-per-cell    ;; {CellX} {Seg} (vectorof Nat) -> {Seg}
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
(define (choose-best-segment-per-column  ;; TM {CellX} {Seg} (vectorof Nat) -> {Seg}
          tm matching-cells all-matching-segments potential-overlaps)
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
(define (get-cells-with-fewest-segments  ;; TM {ColX} -> {CellX}
          tm columnxs)
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
(define (get-active-cells tm)            ;; TM -> {CellX}
  (vector->list (tm-active-cells tm)))
                                                                                            ;
(define (get-predicted-active-cells tm)  ;; TM -> {CellX}
  (tm-predicted-active-cells tm))
                                                                                            ;
(define (get-winner-cells tm)            ;; TM -> {CellX}
  (tm-winner-cells tm))
                                                                                            ;
(define (get-predicted-cells tm)         ;; TM -> {CellX}
  (tm-predicted-cells tm))
                                                                                            ;
(define (get-active-basal-segments tm)   ;; TM -> {CellX}
  (tm-active-basal-segments tm))
                                                                                            ;
(define (get-active-apical-segments tm)  ;; TM -> {CellX}
  (tm-active-apical-segments tm))

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
      (when (fx>=? (syn-perm synapse) threshold)
        (let ((segment    (vector-ref segments segx))
              (new-nacsfs (add1 (vector-ref nacsfs segx))))
          (if (null? reduced-threshold-cells)
              (when (fx=? new-nacsfs actthresh)
                (unless (memq segment actsegs)
                  (set! actsegs (cons segment actsegs))))
              (when (or (and (fx=? new-nacsfs reduced-threshold)
                             (memv (seg-cellx segment) reduced-threshold-cells))
                        (fx=? new-nacsfs actthresh))
                (unless (memq segment actsegs)
                  (set! actsegs (cons segment actsegs)))))
          (vector-set! nacsfs segx new-nacsfs))))
    (define (update-potsegs segx syn-low syn-high)
      (let* ( (synapses (seg-synapses (vector-ref segments segx)))
              (synapse  (search synapses syn-low syn-high)))
        (when synapse
          (let ((new-napsfs (add1 (vector-ref napsfs segx))))
            (when (fx=? new-napsfs potthresh)
              (unless (memq (vector-ref segments segx) potsegs)
                (set! potsegs (cons (vector-ref segments segx) potsegs))))
            (vector-set! napsfs segx new-napsfs))
          (update-actsegs segx synapse))))
    (vector-for-each
      (lambda (cellx)
        (let* ((syn-low  (make-syn cellx min-perm))
               (syn-high (fx+ syn-low max-perm)))
          (for-each
            (lambda (segx)
              (update-potsegs segx syn-low syn-high))
            (vector-ref pre-index cellx))))
      active-presynaptic-cells)
    (values actsegs potsegs napsfs)))
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
(define (destroy-segment tm segment)     ;; TM Segment ->
  ;; make segment's flatx available for re-use
  (let* ( (flatx (seg-flatx segment))
          (null-segment (make-seg -1 flatx (make-synapses 0))))
    (vector-set! (tm-seg-index tm) flatx null-segment)
    (tm-free-flatx-set! tm (cons flatx (tm-free-flatx tm)))))
                                                                                            ;
(define (cellx->colx tm cellx)           ;; TM CellX -> ColX
  (fxdiv cellx (tm-cells-per-column tm)))
                                                                                            ;
(define (add-to-pre-index                ;; CellX Segment (CellVecOf {FlatX}) ->
          prex seg pre-index)
  ;; add to the list of segment flatxs with synapses from the pre-synaptic cell
  (let ((pre-list (vector-ref pre-index prex)))
    (when (or (null? pre-list) (not (fx=? (seg-flatx seg) (car pre-list))))
      (vector-set! pre-index prex (cons (seg-flatx seg) pre-list)))))
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
(define (map-segments-to-cells segments) ;; {Segment} -> {CellX}
  (list-sort fx<? (map seg-cellx segments)))
                                                                                            ;
(define (filter-segments-by-cell         ;; {Seg} {CellX} -> {Seg}
          segments cellxs)
  ;; Return the subset of segments that are on the provided cells.
  (filter
    (lambda (segment)
      (memv (seg-cellx segment) cellxs))
    segments))
                                                                                            ;
(define (map-segments-to-synapse-counts  ;; {Seg} -> {Nat}
          segments)
  (map
    (lambda (segment)
      (synapses-length (seg-synapses segment)))
    segments))
                                                                                            ;
(define (/cpc tm)                        ;; TM -> (CellX -> ColX)
  (lambda (cellx) 
    (fxdiv cellx (tm-cells-per-column tm))))
                                                                                            ;
(define (set-compare tm a b)             ;; TM {CellX} {ColX} -> {CellX} {ColX}
  ;; produce intersection (a&b), difference (b-a), using (a/cpc) for compare
  (let ((a-within-b-mask (in1d (map (/cpc tm) a) b))
        (b-within-a-mask (in1d b (map (/cpc tm) a))))
    (values (include-by-mask a a-within-b-mask) (exclude-by-mask b b-within-a-mask))))
                                                                                           ;
(define (argmax-multi a group-keys)      ;; {Number} {Number} -> {Nat}
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
(define (get-all-cells-in-columns tm     ;; TM {ColX} -> {CellX}
           colxs)
  ;; Calculate all cell indices in the specified columns.
  (let ((cpc (tm-cells-per-column tm)))
    (apply append 
      (map
        (lambda (colx)
          (let ((first (fx* colx cpc)))
            (build-list cpc (lambda (i) (fx+ first i)))))
        colxs))))
                                                                                            ;
(define (prune-synapses syns omits)      ;; Synapses (listof Nat) -> Synapses
  ;; omit from synapses elements indexed by omits (which is sorted)
  (let* ( (reslen (fx- (synapses-length syns) (length omits)))
          (result (make-synapses reslen)))
    (let loop ((rx 0) (sx 0) (omits omits))
      (cond [ (fx=? rx reslen) result ]
            [ (if (null? omits) #f
                (fx=? sx (car omits))) (loop rx (add1 sx) (cdr omits)) ]
            [ else
                (synapses-set! result rx (synapses-ref syns sx))
                (loop (add1 rx) (add1 sx) omits) ]))))
                                                                                            ;
(define  adapt-segment                   ;; TM Seg (vectorof CellX) Connections [(Perm -> Perm) Perm] ->
  ;; Updates synapses on segment. Strengthens active synapses; weakens inactive synapses.
  ;; Remove synapse on zero permanence, destroy segment if no synapses left.
  (case-lambda
  [(tm segment active-input connections)
    ;; use increment/decrement parameter values
    (adapt-segment tm segment active-input connections
                   (lambda (p) (clip-max (fx+ p (tm-permanence-increment tm))))
                   (tm-permanence-decrement tm))]
  [(tm segment active-input connections increment-proc permanence-decrement)
    ;; general case: increment supplied by proc
    (let ((synapses (seg-synapses segment)))
      (let build-s-t-d ((sx (fx- (synapses-length synapses) 1))
                        (synapses-to-destroy '()))
        (if (negative? sx)
            (cond [ (null? synapses-to-destroy) ]
                  [ (fx=? (length synapses-to-destroy) (synapses-length synapses))
                      (vector-set! connections (seg-cellx segment)
                        (remq segment (vector-ref connections (seg-cellx segment))))
                      (destroy-segment tm segment) ]
                  [ else (seg-synapses-set! segment 
                           (prune-synapses synapses synapses-to-destroy)) ])
            (let* ( (synapse    (synapses-ref synapses sx))
                    (prex       (syn-prex synapse))
                    (permanence (if (search active-input prex prex)
                                  (increment-proc (syn-perm synapse))
                                  (clip-min (fx- (syn-perm synapse) permanence-decrement)))))
              (if (zero? permanence)
                  ;; build synapses-to-destroy indices as sorted list
                  (build-s-t-d (fx- sx 1) (cons sx synapses-to-destroy))
                  (begin
                    (synapses-set! synapses sx (make-syn prex permanence))
                    (build-s-t-d (fx- sx 1) synapses-to-destroy)))))))]))
                                                                                            ;
(define (adjust-synapses tm              ;; TM Connections {Seg} (vectorof CellX) ->
          connections segments active-input)
  ;; For each specified segment, update the permanence of each synapse
  ;; according to whether the synapse would be active given the specified
  ;; active inputs.
  (for-each
    (lambda (segment)
      (adapt-segment tm segment (list->vector active-input) connections))
    segments))
                                                                                            ;
(define (adjust-active-synapses tm       ;; TM Connections {Seg} {CellX} Permanence ->
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
(define (grow-synapses-to-sample tm      ;; TM {Seg} {CellX} Nat|{Nat} (CellVecOf {FlatX}) ->
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
(define (in-synapses? cellx synapses)    ;; CellX (vectorof Synapse) -> Boolean
  ;; produce whether cellx equals any pre-synaptic-cell index of synapses
  (let* ( (syn-low  (make-syn cellx min-perm))
          (syn-high (fx+ syn-low max-perm)))
    (search synapses syn-low syn-high)))
                                                                                            ;
(define (grow-synapses                   ;; TM Seg Nat (listof CellX) (CellVecOf {FlatX}) ->
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
        (let ((new-prexs (vector-sample (list->vector candidates) n-actual))
              (init-perm (tm-initial-permanence tm)))
          (seg-synapses-set! segment     ;; append new synapses to segment
            (build-synapses (fx+ num-synapses n-actual)
              (lambda (sx)
                (if (fx<? sx num-synapses)
                    (synapses-ref synapses sx)
                    (let ((new-prex (vector-ref new-prexs (fx- sx num-synapses))))
                      (add-to-pre-index new-prex segment pre-index)
                      (make-syn new-prex init-perm))))))
          (vector-sort! fx<? (seg-synapses segment))))))
                                                                                            ;
(define (create-segment tm               ;; TM CellX Connections -> Seg
          cellx connections)
  ;; produce a new segment on the cell, updating the index of segments
  (let* ( (segments (vector-ref connections cellx))
          (new-flatx (if (null? (tm-free-flatx tm))
                          (let ((next (tm-next-flatx tm)))
                            (when (fx=? next (vector-length (tm-seg-index tm)))
                              (tm-seg-index-set! tm (vector-extend (tm-seg-index tm))))
                            next)  
                          (car (tm-free-flatx tm))))
          (segment (make-seg cellx new-flatx (make-synapses 0))))
    (vector-set! (tm-seg-index tm) new-flatx segment)
    (if (null? (tm-free-flatx tm))
      (tm-next-flatx-set! tm (add1 new-flatx))
      (tm-free-flatx-set! tm (cdr (tm-free-flatx tm))))
    (vector-set! connections cellx (cons segment segments))
    segment))
                                                                                            ;
(define (create-segments tm              ;; TM Connections {CellX} -> {Seg}
          connections cellxs)
  (map
    (lambda (cellx)
      (create-segment tm cellx connections))
    cellxs))
                                                                                            ;
)