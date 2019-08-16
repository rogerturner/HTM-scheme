#!chezscheme

;; === HTM-scheme Apical Tiebreak Temporal Memory algorithm Copyright 2019 Roger Turner. ===
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
  #|

A generalized Temporal Memory with apical dendrites that add a "tiebreak".

Basal connections are used to implement traditional Temporal Memory.

The apical connections are used for further disambiguation. If multiple cells
in a minicolumn have active basal segments, each of those cells is predicted,
unless one of them also has an active apical segment, in which case only the
cells with active basal and apical segments are predicted.

In other words, the apical connections have no effect unless the basal input
is a union of SDRs (e.g. from bursting minicolumns).
  |#

(library (HTM-scheme HTM-scheme algorithms apical_tiebreak_temporal_memory)
                                                                                            ;
(export
  reset
  depolarize-cells
  activate-cells
  (rename
    (tm-active-cells           get-active-cells)
    (tm-predicted-active-cells get-predicted-active-cells)
    (tm-winner-cells           get-winner-cells)
    (tm-predicted-cells        get-predicted-cells)
    (tm-active-basal-segments  get-active-basal-segments)
    (tm-active-apical-segments get-active-apical-segments)
    (tm-n-segments-created     get-n-segments-created)
    (tm-n-synapses-created     get-n-synapses-created))
  tm                                     ;; for pair/sequence memory
  make-tm
  number-of-cells
  map-segments-to-cells
  cols-from-cells                        ;; for l2l4_patch
  (rename                                ;; for temporal_memory_test
    (create-segment       test:create-segment)
    (tm-basal-connections test:tm-basal-connections)
    (grow-synapses        test:grow-synapses)
    (tm-basal-pre-index   test:tm-basal-pre-index)
    (tm-next-flatx        test:tm-next-flatx)))
                                                                                            ;
(import
  (except (chezscheme) add1 make-list random reset)
  (HTM-scheme HTM-scheme algorithms htm_prelude)
  (HTM-scheme HTM-scheme algorithms htm_concept))
  
;; === Types (see htm_concept.ss) ===
                                                                                            ;
;; FlatX        = Nat, segment sequence number
;; CellX        = Nat, cell index
;; CellVecOf    = Vector indexed by CellX
;; ColX         = Nat, column index of cell
;; TM           = Record: tm parameters, cells

;; === Layer record ===
                                                                                            ;
(define-record-type tm                   ;; TM
  (fields
    column-count                         ;; Nat
    basal-input-size                     ;; Nat
    apical-input-size                    ;; Nat
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
    (mutable active-cells)               ;; (listof CellX)
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
    (mutable free-flatx)                 ;; (listof Nat): indices available for re-use
    (mutable n-segments-created)
    (mutable n-synapses-created))
  (protocol
    (lambda (new)                        ;; (listof KWarg) -> TM
      (lambda (kwargs)
        (let* ( (tm (apply new (key-word-args kwargs attm-defaults)))
                (num-cells (* (tm-column-count tm) (tm-cells-per-column tm))))
          (tm-basal-connections-set!  tm (make-vector num-cells '()))
          (tm-apical-connections-set! tm (make-vector num-cells '()))
          (tm-basal-pre-index-set!    tm (make-vector num-cells '()))
          (tm-apical-pre-index-set!   tm (make-vector (tm-apical-input-size tm) '()))
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
  [active-cells                          . ()]
  [winner-cells                          . ()]
  [predicted-cells                       . ()]
  [predicted-active-cells                . ()]
  [active-basal-segments                 . ()]
  [active-apical-segments                . ()]
  [matching-basal-segments               . ()]
  [matching-apical-segments              . ()]
  [basal-potential-overlaps              . #vfx()]
  [apical-potential-overlaps             . #vfx()]
  [basal-pre-index                       . #()]
  [apical-pre-index                      . #()]
  [seg-index                             . #()]
  [next-flatx                            . 0]
  [free-flatx                            . ()]
  [n-segments-created                    . 0]
  [n-synapses-created                    . 0] ))

;; === Apical Tiebreak Temporal Memory Algorithm ===
                                                                                            ;
(define (reset tm)                       ;; TM ->
  ;; clear all cell and segment activity
  (tm-active-cells-set!              tm '())
  (tm-winner-cells-set!              tm '())
  (tm-predicted-cells-set!           tm '())
  (tm-predicted-active-cells-set!    tm '())
  (tm-active-basal-segments-set!     tm '())
  (tm-active-apical-segments-set!    tm '())
  (tm-matching-basal-segments-set!   tm '())
  (tm-matching-apical-segments-set!  tm '())
  (tm-basal-potential-overlaps-set!  tm (fxvector))
  (tm-apical-potential-overlaps-set! tm (fxvector)))
                                                                                            ;
(define (depolarize-cells                ;; TM {CellX} {CellX} Boolean ->
          tm basal-input apical-input learn)
  ;; calculate predictions, save active/matching segments
  (let-values ([(active-apical-segments matching-apical-segments apical-potential-overlaps)
                (calculate-apical-segment-activity tm apical-input)])
    (let ((reduced-basal-threshold-cells
            (if (or learn (not (tm-use-apical-modulation-basal-threshold tm)))
              '()
              (map-segments-to-cells active-apical-segments))))
      (let-values ([(active-basal-segments matching-basal-segments basal-potential-overlaps)
                    (calculate-basal-segment-activity tm basal-input
                        reduced-basal-threshold-cells)])
        (tm-predicted-cells-set! tm 
          (calculate-predicted-cells tm active-basal-segments active-apical-segments))
        (tm-active-basal-segments-set!     tm active-basal-segments)
        (tm-active-apical-segments-set!    tm active-apical-segments)
        (tm-matching-basal-segments-set!   tm matching-basal-segments)
        (tm-matching-apical-segments-set!  tm matching-apical-segments)
        (tm-basal-potential-overlaps-set!  tm basal-potential-overlaps)
        (tm-apical-potential-overlaps-set! tm apical-potential-overlaps)))))
                                                                                            ;
(define (activate-cells                  ;; TM {ColX} {CellX} {CellX} {CellX} {CellX} Boolean ->
          tm feedforward-input basal-reinforce-candidates apical-reinforce-candidates
          basal-growth-candidates apical-growth-candidates learning . bursting-columns)
  ;; activate cells in cols with feedforward-input; cols with no predicted cell burst; learn
  (let* (                                ;; calculate active cells
      (correct-predicted-cells 
        (cells-in-cols tm (tm-predicted-cells tm) feedforward-input))
      (bursting-columns   
        (if (pair? bursting-columns)
          (car bursting-columns)
          (setdiff1d feedforward-input (cols-from-cells tm (tm-predicted-cells tm)))))
      (new-active-cells 
        (append correct-predicted-cells (get-all-cells-in-columns tm bursting-columns))))
  (let*-values (                       ;; calculate learning
        [ ( learning-active-basal-segments
            learning-matching-basal-segments
            basal-segments-to-punish
            new-basal-segment-cells
            learning-cells)
              (calculate-basal-learning tm
                feedforward-input bursting-columns correct-predicted-cells
                (tm-active-basal-segments tm) (tm-matching-basal-segments tm)
                (tm-basal-potential-overlaps tm)) ]
        [ ( learning-active-apical-segments
            learning-matching-apical-segments
            apical-segments-to-punish
            new-apical-segment-cells)
              (calculate-apical-learning tm
                learning-cells feedforward-input
                (tm-active-apical-segments tm) (tm-matching-apical-segments tm)
                (tm-apical-potential-overlaps tm)) ] )
      (when learning
        ;; Learn on existing segments
        (let ((basal-reinforce-candidates  (list->vector basal-reinforce-candidates))
              (apical-reinforce-candidates (list->vector apical-reinforce-candidates)))
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
                                    (- (tm-apical-predicted-segment-decrement tm)))))
        ;; Grow new segments
        (when (positive? (length basal-growth-candidates))
          (learn-on-new-segments tm (tm-basal-connections tm) new-basal-segment-cells 
                                 basal-growth-candidates (tm-basal-pre-index tm)))
        (when (positive? (length apical-growth-candidates))
          (learn-on-new-segments tm (tm-apical-connections tm) new-apical-segment-cells
                                 apical-growth-candidates (tm-apical-pre-index tm))))
      ;; Save the results
      (tm-active-cells-set! tm (sort! fx<? new-active-cells))
      (tm-winner-cells-set! tm learning-cells)
      (tm-predicted-active-cells-set! tm correct-predicted-cells))))
                                                                                          ;
(define (calculate-basal-learning tm     ;; TM {ColX} {ColX} {CellX} {Seg} {Seg} {Nat}
          feedforward-input bursting-columns correct-predicted-cells
          active-basal-segments matching-basal-segments basal-potential-overlaps)
  ;; Basic Temporal Memory learning. Correctly predicted cells always have
  ;; active basal segments, and we learn on these segments. In bursting
  ;; columns, we either learn on an existing basal segment, or we grow a new one.
  ;; The only influence apical dendrites have on basal learning is: the apical
  ;; dendrites influence which cells are considered "predicted". So an active
  ;; apical dendrite can prevent some basal segments in active columns from
  ;; learning.
  (let* ( 
      (learning-active-basal-segments 
        (filter-segments-by-cell active-basal-segments correct-predicted-cells))
      (cells-for-matching-basal
        (map-segments-to-cells matching-basal-segments))
      (basal-segments-to-punish      ;; here because unique! will mutate cells-for-matching-basal
        (exclude-segments tm matching-basal-segments cells-for-matching-basal feedforward-input))
      (matching-cells (unique! fx=? cells-for-matching-basal))
      (matching-cells-in-bursting-columns
        (cells-in-cols tm matching-cells bursting-columns))
      (bursting-columns-with-no-match
        (setdiff1d bursting-columns (cols-from-cells tm matching-cells)))
      (learning-matching-basal-segments
        (choose-best-segment-per-column tm matching-cells-in-bursting-columns
          matching-basal-segments basal-potential-overlaps))
      (new-basal-segment-cells
        (get-cells-with-fewest-segments tm bursting-columns-with-no-match))
      (learning-cells
        (unique! fx=? (sort! fx<?
          (append correct-predicted-cells
                  (map-segments-to-cells learning-matching-basal-segments)
                  new-basal-segment-cells 
                  '())))))       ;; *copy new-basal-segment-cells*
    (values
      learning-active-basal-segments
      learning-matching-basal-segments
      basal-segments-to-punish
      new-basal-segment-cells
      learning-cells)))
                                                                                          ;
(define (calculate-apical-learning tm    ;; TM {CellX} {ColX} {Seg} {Seg} {Nat}
          learning-cells feedforward-input  
          active-apical-segments matching-apical-segments apical-potential-overlaps)
  ;; Calculate apical learning for each learning cell.
  ;; The set of learning cells was determined completely from basal segments.
  ;; Do all apical learning on the same cells.
  ;; Learn on any active segments on learning cells. For cells without active
  ;; segments, learn on the best matching segment. For cells without a matching
  ;; segment, grow a new segment.
  (let* ( 
      (learning-active-apical-segments 
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
      (apical-segments-to-punish
        (exclude-segments tm matching-apical-segments cells-for-matching-apical feedforward-input)))
    (values
      learning-active-apical-segments
      learning-matching-apical-segments
      apical-segments-to-punish
      new-apical-segment-cells)))
                                                                                          ;
(define (calculate-apical-segment-activity ;; TM {CellX} -> {Seg} {Seg} (vectorof Nat)
          tm active-input)
  ;; calculate active/matching apical segments
  (calculate-segment-activity tm active-input (tm-apical-pre-index tm) 0 '()))
                                                                                            ;
(define (calculate-basal-segment-activity  ;; TM {CellX} {CellX} -> {Seg} {Seg} (vectorof Nat)
          tm active-input reduced-basal-threshold-cells)
  ;; calculate active/matching basal segments
  (calculate-segment-activity tm active-input (tm-basal-pre-index tm)
    (tm-reduced-basal-threshold tm) reduced-basal-threshold-cells))
                                                                                            ;
(define (calculate-predicted-cells tm    ;; TM {Seg} {Seg} -> {CellX}
          active-basal-segments active-apical-segments)
  ;; calculate predicted cells for given active segments
  (let ((cells-for-basal-segments (unique! fx=? (map-segments-to-cells active-basal-segments))))
    (if (tm-use-apical-tiebreak tm)
      (let ((cells-for-apical-segments (unique! fx=? (map-segments-to-cells active-apical-segments))))
        (let loop ( (cfbs  cells-for-basal-segments )
                    (cfas  cells-for-apical-segments)
                    (fully-depolarized-cells  (list))
                    (partly-depolarized-cells (list))
                    (inhibited-cols           (list)))
          (if (null? cfbs)
            (sort! fx<?
              (append!
                (filter (lambda (cellx)
                    (not (memv (cellx->colx tm cellx) inhibited-cols)))
                  partly-depolarized-cells)
                fully-depolarized-cells))
            (let ((cell-with-basal-seg (car cfbs))
                  (cell-with-apical-seg 
                    (if (null? cfas)  (greatest-fixnum)
                      (car cfas))))
              (cond
                [(fx=? cell-with-basal-seg cell-with-apical-seg)   ;; fully depolarized cell
                  (loop (cdr cfbs) (cdr cfas)
                    (cons cell-with-basal-seg fully-depolarized-cells)
                    partly-depolarized-cells
                    (cons (cellx->colx tm cell-with-basal-seg) inhibited-cols)) ]
                [(fx<? cell-with-basal-seg cell-with-apical-seg)   ;; partly depolarized cell
                  (loop (cdr cfbs) cfas
                    fully-depolarized-cells
                    (cons cell-with-basal-seg partly-depolarized-cells)
                    inhibited-cols) ]
                [ else                                             ;; apical only
                  (loop cfbs (cdr cfas) fully-depolarized-cells partly-depolarized-cells inhibited-cols)])))))
      cells-for-basal-segments)))
                                                                                            ;
(define (learn tm connections            ;; TM Connections {Seg} {CellX} {CellX} (vectorof Nat) (CellVecOf {FlatX}) ->
          learning-segments active-input growth-candidates potential-overlaps pre-index)
  ;; adjust synapse permanences, grow new synapses, and grow new segments
  (adjust-synapses tm connections learning-segments active-input)
  (let* ( 
      (ss (tm-sample-size tm))
      (max-new                           ;; Nat | {Nat}
        (if (fxnegative? ss)
          (length growth-candidates)
          (map (lambda (seg)
              (fx- ss (fxvector-ref potential-overlaps (seg-flatx seg))))
            learning-segments)))
      (msps (tm-max-synapses-per-segment tm))
      (max-new
          (if (fxnegative? msps)
            max-new
            (if (fxnegative? ss)
              (map (lambda (seg)
                  (fxmin max-new
                    (fx- msps (synapses-length (seg-synapses seg)))))
                learning-segments)
              (map (lambda (seg m-n)
                  (fxmin m-n
                    (fx- msps (synapses-length (seg-synapses seg)))))
                learning-segments max-new)))))
    (grow-synapses-to-sample tm learning-segments growth-candidates max-new pre-index)))
                                                                                            ;
(define (learn-on-new-segments tm        ;; TM Connections {CellX} {CellX} (CellVecOf {FlatX}) ->
          connections new-segment-cells growth-candidates pre-index)
  ;; create segments for new-segment-cells and connect to growth-candidates
  (let* ( 
      (num-new-synapses (length growth-candidates))
      (num-new-synapses (if (fxnegative? (tm-sample-size tm))
                          num-new-synapses
                          (fxmin num-new-synapses (tm-sample-size tm))))
      (num-new-synapses (if (fxnegative? (tm-max-synapses-per-segment tm))
                          num-new-synapses
                          (fxmin num-new-synapses (tm-max-synapses-per-segment tm))))
      (new-segments (create-segments tm connections new-segment-cells)))
    (grow-synapses-to-sample tm new-segments growth-candidates num-new-synapses pre-index)))
                                                                                            ;
(define (choose-best-segment-per-cell    ;; {CellX} {Seg} (vectorof Nat) -> {Seg}
          cells all-matching-segments potential-overlaps)
  ;; choose matching segment with max active potential synapses (on tie choose first)
  (let next-cell ((cells cells) (learning-segments '()))
    (if (null? cells) 
      (reverse learning-segments)
      (let next-segment ( (segments all-matching-segments) 
                          (best-seg #f) 
                          (max-overlap (least-fixnum)))
        (cond
          [ (null? segments) (next-cell (cdr cells) (cons best-seg learning-segments)) ]
          [ (fx=? (car cells) (seg-cellx (car segments)))
              (let ((overlap (fxvector-ref potential-overlaps (seg-flatx (car segments)))))
                (if (fx>? overlap max-overlap)
                    (next-segment (cdr segments) (car segments) overlap)
                    (next-segment (cdr segments) best-seg max-overlap))) ]
          [ else (next-segment (cdr segments) best-seg max-overlap) ] )))))
                                                                                            ;
(define (choose-best-segment-per-column  ;; TM {CellX} {Seg} (vectorof Nat) -> {Seg}
          tm matching-cells all-matching-segments potential-overlaps)
  ;; like choose-best-segment-per-cell but for cells in col
  (let* ( 
      (candidate-segments
        (filter-segments-by-cell all-matching-segments matching-cells))
      (cell-scores
        (map (lambda (segment)
            (fxvector-ref potential-overlaps (seg-flatx segment)))
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
  ;; produce cell for each col with fewest total basal segments; break ties randomly
  (map (lambda (colx)
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
(define (number-of-cells tm)             ;; TM -> Nat
  (fx* (tm-column-count tm) (tm-cells-per-column tm)))

;; === Supporting Functions ===
                                                                                            ;
;; --- Connections: see NuPIC connections.py ---
                                                                                            ;
(define (calculate-segment-activity      ;; TM {CellX} (CellVecOf {FlatX}) Nat {CellX} -> {Segment} {Segment} (Vectorof Nat)
          tm active-input pre-index reduced-threshold reduced-threshold-cells)
  ;; for active-input cells, produce lists of segments with connected/potential synapses,
  ;; and potential synapse counts (combines connections.computeActivity and _calculateApical/BasalSegmentActivity)
  (let (
    (napsfs     (make-fxvector (tm-next-flatx tm) 0)) ;; "num-active-potential-synapses-for-segment"
    (nacsfs     (make-fxvector (tm-next-flatx tm) 0)) ;; "num-active-connected-synapses-for-segment"
    (threshold  (tm-connected-permanence tm))
    (segments   (tm-seg-index tm))
    (actsegs    '())
    (potsegs    '())
    (actthresh  (tm-activation-threshold tm))
    (potthresh  (tm-min-threshold tm))
    (no-reduced (null? reduced-threshold-cells)))
    (define (update-actsegs segx synapse)
      (when (fx>=? (syn-perm synapse) threshold)
        (let ((segment    (vector-ref segments segx))
              (new-nacsfs (add1 (fxvector-ref nacsfs segx))))
          (if no-reduced
              (when (fx=? new-nacsfs actthresh)
                (unless (memq segment actsegs)
                  (set! actsegs (cons segment actsegs))))
              (when (or (and (fx=? new-nacsfs reduced-threshold)
                             (memv (seg-cellx segment) reduced-threshold-cells))
                        (fx=? new-nacsfs actthresh))
                (unless (memq segment actsegs)
                  (set! actsegs (cons segment actsegs)))))
          (fxvector-set! nacsfs segx new-nacsfs))))
    (define (update-potsegs segx syn-low syn-high)
      (let* ( (synapses (seg-synapses (vector-ref segments segx)))
              (synapse  (synapses-search synapses syn-low syn-high)))
        (when synapse
          (let ((new-napsfs (add1 (fxvector-ref napsfs segx))))
            (when (fx=? new-napsfs potthresh)
              (set! potsegs (cons (vector-ref segments segx) potsegs)))
            (fxvector-set! napsfs segx new-napsfs))
          (update-actsegs segx synapse))))
    (for-each (lambda (cellx)
        (let* ((syn-low  (make-syn cellx min-perm))
               (syn-high (fx+ syn-low max-perm)))
          (for-each (lambda (segx)
              (update-potsegs segx syn-low syn-high))
            (vector-ref pre-index cellx))))
      active-input)
    (values actsegs
      (unique! eq? (sort!
        (lambda (sega segb)
          (cond [(fx<? (seg-cellx sega) (seg-cellx segb))]
                [(fx=? (seg-cellx sega) (seg-cellx segb))
                    (fx<? (seg-flatx sega) (seg-flatx segb))]
                [else #f]))
        potsegs))
      napsfs)))
                                                                                            ;
(define (map-segments-to-cells segments) ;; {Seg} -> {CellX}
  ;; produce cellxs from segments
  (sort! fx<? (map seg-cellx segments)))
                                                                                            ;
(define (filter-segments-by-cell         ;; {Seg} {CellX} -> {Seg}
          segments cellxs)
  ;; produce subset of segments that are on cellxs
  (filter
    (lambda (segment)
      (memv (seg-cellx segment) cellxs))
    segments))
                                                                                            ;
(define (map-segments-to-synapse-counts  ;; {Seg} -> {Nat}
          segments)
  (map (lambda (segment)
      (synapses-length (seg-synapses segment)))
    segments))
                                                                                            ;
(define (cellx->colx tm cellx)           ;; TM CellX -> ColX
  (fxdiv cellx (tm-cells-per-column tm)))
                                                                                            ;
(define (/cpc tm)                        ;; TM -> (CellX -> ColX)
  (lambda (cellx) (cellx->colx tm cellx)))
                                                                                            ;
(define (cols-from-cells tm cellxs)      ;; TM {CellX} -> {ColX}
  ;; produce canonical colxs of the cellxs (which are sorted)
  (unique! fx=? (map (/cpc tm) cellxs)))
                                                                                            ;
(define (cells-in-cols tm cellxs colxs)  ;; TM {CellX} {ColX} -> {CellX}
  ;; produce cellxs for which (col of) cellx is in colxs
  (reverse!
    (fold-left (lambda (acc cellx)
        (if (memv (cellx->colx tm cellx) colxs)
          (cons cellx acc)
          acc))
      (list)
      cellxs)))
                                                                                            ;
(define (exclude-segments tm             ;; TM {Seg} {CellX} {ColX} -> {Seg}
          segs cellxs colxs)
  ;; produce segments for which cellx is not in colxs; segs and cellxs are parallel
  (fold-left (lambda (acc seg cellx)
      (if (memv (cellx->colx tm cellx) colxs)
        acc
        (cons seg acc)))
    (list)
    segs cellxs))
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
  ;; produce all cellxs in the colxs
  (let ((cpc (tm-cells-per-column tm)))
    (apply append! 
      (map (lambda (colx)
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
(define (destroy-segment tm segment)     ;; TM Seg ->
  ;; make segment's flatx available for re-use
  (let* ( (flatx (seg-flatx segment))
          (null-segment (make-seg -1 flatx (make-synapses 0))))
    (vector-set! (tm-seg-index tm) flatx null-segment)
    (tm-free-flatx-set! tm (cons flatx (tm-free-flatx tm)))))
                                                                                            ;
(define (adapt-segment tm segment        ;; TM Seg (vectorof CellX) Connections Perm Perm ->
          active-input connections permanence-delta permanence-decrement)
  ;; Updates synapses on segment. Strengthens active synapses; weakens inactive synapses.
  ;; Remove synapse on zero permanence, destroy segment if no synapses left.
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
                  (permanence (if (fxsearch active-input prex)
                                (clip-max (clip-min (fx+ (syn-perm synapse) permanence-delta)))
                                (clip-min (fx- (syn-perm synapse) permanence-decrement)))))
            (if (zero? permanence)
                ;; build synapses-to-destroy indices as sorted list
                (build-s-t-d (fx- sx 1) (cons sx synapses-to-destroy))
                (begin
                  (synapses-set! synapses sx (make-syn prex permanence))
                  (build-s-t-d (fx- sx 1) synapses-to-destroy))))))))
                                                                                            ;
(define (adjust-synapses tm              ;; TM Connections {Seg} {CellX} ->
          connections segments active-input)
  ;; For each specified segment, update the permanence of each synapse
  ;; according to whether the synapse would be active given the specified
  ;; active inputs.
  (let ((inc (tm-permanence-increment tm))
        (dec (tm-permanence-decrement tm)))
    (for-each (lambda (segment)
        (adapt-segment tm segment active-input connections inc dec))
      segments)))
                                                                                            ;
(define (adjust-active-synapses tm       ;; TM Connections {Seg} {CellX} Perm ->
          connections segments active-input permanence-delta)
  ;; For each specified segment, add a delta to the permanences of the
  ;; synapses that would be active given the specified active inputs.
  (for-each (lambda (segment)
      (adapt-segment tm segment active-input connections permanence-delta 0))
    segments))
                                                                                          ;
(define (grow-synapses-to-sample tm      ;; TM {Seg} {CellX} Nat|{Nat} (CellVecOf {FlatX}) ->
          segments inputs sample-size pre-index)
  ;; For each segment, grow synapses to a random subset of the
  ;; inputs that aren't already connected to the segment.
  (if (fixnum? sample-size)
    (for-each (lambda (segment)
        (grow-synapses tm segment sample-size inputs pre-index))
      segments)
    (for-each (lambda (segment sample-size)
        (grow-synapses tm segment sample-size inputs pre-index))
      segments sample-size)))
                                                                                            ;
(define (in-synapses? cellx synapses)    ;; CellX Synapses -> Boolean
  ;; produce whether cellx equals any pre-synaptic-cell index of synapses
  (let* ( (syn-low  (make-syn cellx min-perm))
          (syn-high (fx+ syn-low max-perm)))
    (synapses-search synapses syn-low syn-high)))
                                                                                            ;
(define (add-to-pre-index                ;; CellX Seg (CellVecOf {FlatX}) ->
          prex seg pre-index)
  ;; add to the list of segment flatxs with synapses from the pre-synaptic cell
  (let ((pre-list (vector-ref pre-index prex)))
    (when (or (null? pre-list) (not (fx=? (seg-flatx seg) (car pre-list))))
      (vector-set! pre-index prex (cons (seg-flatx seg) pre-list)))))
                                                                                            ;
(define (grow-synapses                   ;; TM Seg Nat {CellX} (CellVecOf {FlatX}) ->
          tm segment n-desired-new-synapses inputs pre-index)
  ;; grow synapses to all inputs that aren't already connected to the segment
  (let build-candidates ((cs '()) (is inputs) (nc 0))
    (cond 
      [ (null? is)
        (if (fx<=? nc n-desired-new-synapses)
          (let ((init-perm (tm-initial-permanence tm)))
            (tm-n-synapses-created-set! tm (fx+ (tm-n-synapses-created tm) nc))
            (seg-synapses-set! segment 
              (list->synapses (merge! fx<? (sort! fx<?
                                              (let map-car! ((c cs))
                                                (cond [(null? c) cs]
                                                  [else
                                                    (add-to-pre-index (car c) segment pre-index)
                                                    (set-car! c (make-syn (car c) init-perm))
                                                    (map-car! (cdr c)) ])))
                                           (synapses->list (seg-synapses segment))))))
          (when (positive? n-desired-new-synapses)
            (tm-n-synapses-created-set! tm (fx+ (tm-n-synapses-created tm) n-desired-new-synapses))
            (let ((new-prexs (vector-sample (list->vector cs) n-desired-new-synapses))
                  (init-perm (tm-initial-permanence tm)))
              (vector-sort! fx<? new-prexs)
              (seg-synapses-set! segment 
                (list->synapses (merge! fx<? (build-list n-desired-new-synapses (lambda (newx) 
                                                (let ((new-prex (vector-ref new-prexs newx)))
                                                  (add-to-pre-index new-prex segment pre-index)
                                                  (make-syn new-prex init-perm))))
                                             (synapses->list (seg-synapses segment)))))))) ]
      [ (in-synapses? (car is) (seg-synapses segment))
          (build-candidates cs (cdr is) nc) ]
      [ else (build-candidates (cons (car is) cs) (cdr is) (add1 nc)) ])))
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
    (tm-n-segments-created-set! tm (add1 (tm-n-segments-created tm)))
    segment))
                                                                                            ;
(define (create-segments tm              ;; TM Connections {CellX} -> {Seg}
          connections cellxs)
  (map (lambda (cellx)
      (create-segment tm cellx connections))
    cellxs))
                                                                                            ;
)