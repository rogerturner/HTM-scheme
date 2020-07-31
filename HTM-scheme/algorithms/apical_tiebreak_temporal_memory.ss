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
  #|

Translated from numenta htmresearch/.../apical_tiebreak_temporal_memory.py,
htmresearch/.../numpy_helpers.py, nupic-core/.../SparseMatrixConnections.cpp --
see comments there for descriptions of functions and parameters.
Follows logic of apical_tiebreak_temporal_memory.py except:
  active-apical-segments always used to produce reduced-basal-threshold-cells
  activate-cells has an optional bursting-columns input
  separate growth-candidates inputs for cortical column and minicolumn
(see layer4.ss for usage of latter 2 features)

See htm_concept.ss for type and data structure description and code conventions.
  
Description from apical_tiebreak_temporal_memory.py:
  A generalized Temporal Memory with apical dendrites that add a "tiebreak".
  Basal connections are used to implement traditional Temporal Memory.
  The apical connections are used for further disambiguation. If multiple cells
  in a minicolumn have active basal segments, each of those cells is predicted,
  unless one of them also has an active apical segment, in which case only the
  cells with active basal and apical segments are predicted.

    |#

(library (HTM-scheme HTM-scheme algorithms apical_tiebreak_temporal_memory)
                                                                                            ;
(export
make-tm
reset
depolarize-cells
activate-cells
  (rename
    (tm-active-cells           get-active-cells)
    (tm-predicted-active-cells get-predicted-active-cells)
    (tm-winner-cells           get-winner-cells)
    (tm-predicted-cells        get-predicted-cells)
    (tm-active-basal-segments  get-active-basal-segments)
    (tm-active-apical-segments get-active-apical-segments))
  tm                                     ;; for pair/sequence memory
  number-of-cells
  number-of-basal-segments
  number-of-basal-synapses
  number-of-apical-segments
  number-of-apical-synapses
  average-preindex-length
  map-segments-to-cells
  cols-from-cells                        ;; for l2l4_patch
  get-all-cells-in-columns
  (rename                                ;; for temporal_memory_test
    (create-segment       test:create-segment)
    (tm-basal-connections test:tm-basal-connections)
    (grow-synapses        test:grow-synapses)
    (tm-basal-pre-index   test:tm-basal-pre-index)
    (tm-apical-pre-index  test:tm-apical-pre-index)))
                                                                                            ;
(import
  (except (chezscheme) reset)
  (HTM-scheme HTM-scheme algorithms htm_prelude)
  (HTM-scheme HTM-scheme algorithms htm_concept))

;; === Types (see htm_concept.ss) ===
                                                                                            ;
;; CellX        = Nat, cell index
;; CellVecOf    = Vector indexed by CellX
;; ColX         = Nat, column index of cell
;; PreIndex     = Vectorof {Segment}
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
    use-apical-modulation-basal-threshold;; Boolean
    perm-trim-threshold                  ;; Perm
    has-recurrent-connections            ;; Boolean
    same-col?                            ;; (lambda (CellX InputX))
    (mutable apical-activation-threshold);; Nat
    (mutable apical-sample-size)         ;; Nat
    (mutable iteration)                  ;; Fixnum [incrementing by #x10000]
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
    (mutable basal-pre-index)            ;; (vectorof {Segment})
    (mutable apical-pre-index))          ;; (vectorof {Segment})
  (protocol
    (lambda (new)                        ;; (listof KWarg) -> TM
      (lambda (kwargs)
        (let* ( [tm (apply new (key-word-args kwargs attm-defaults))]
                [num-cells (* (tm-column-count tm) (tm-cells-per-column tm))])
          (unless (tm-apical-activation-threshold tm)
            (tm-apical-activation-threshold-set! tm (tm-activation-threshold tm)))
          (unless (tm-apical-sample-size tm)
            (tm-apical-sample-size-set! tm (tm-sample-size tm)))
          (tm-basal-connections-set!  tm (make-vector num-cells '()))
          (tm-apical-connections-set! tm (make-vector num-cells '()))
          (tm-basal-pre-index-set!    tm (make-vector (tm-basal-input-size  tm) '()))
          (tm-apical-pre-index-set!   tm (make-vector (tm-apical-input-size tm) '()))
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
  [perm-trim-threshold                   . ,(perm 0.0)]
  [has-recurrent-connections             . #f]
  [same-col?                             . ,(lambda _ #t)]
  [apical-activation-threshold           . #f]
  [apical-sample-size                    . #f]
  [iteration                             . 0]
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
  [basal-pre-index                       . #()]
  [apical-pre-index                      . #()] ))

;; === Apical Tiebreak Temporal Memory Algorithm ===
                                                                                            ;
(define (reset tm)                       ;; TM ->
  ;; clear all cell and segment activity (overlap counts are stored in segments)
  (tm-active-cells-set!              tm '())
  (tm-winner-cells-set!              tm '())
  (tm-predicted-cells-set!           tm '())
  (tm-predicted-active-cells-set!    tm '())
  (tm-active-basal-segments-set!     tm '())
  (tm-active-apical-segments-set!    tm '())
  (tm-matching-basal-segments-set!   tm '())
  (tm-matching-apical-segments-set!  tm '()))
                                                                                            ;
(define (depolarize-cells                ;; TM {CellX} {CellX} Boolean ->
          tm basal-input apical-input learn)
  ;; calculate predictions, save predicted cells and active/matching segments
  (tm-iteration-set! tm (fx+ #x10000 (tm-iteration tm)))
  (let-values (
      [(active-apical-segments matching-apical-segments)
        (calculate-segment-activity tm
          apical-input (tm-apical-pre-index tm) '() (tm-apical-activation-threshold tm))])
    (let ([reduced-basal-threshold-cells ;; **htmresearch/attm.py has () if learn**
            (if (tm-use-apical-modulation-basal-threshold tm)
              (unique! fx=? (map-segments-to-cells active-apical-segments))
              '())])
      (let-values (
          [(active-basal-segments matching-basal-segments)
            (calculate-segment-activity tm
              basal-input (tm-basal-pre-index tm) reduced-basal-threshold-cells
              (tm-activation-threshold tm))])
        (tm-predicted-cells-set! tm
          (calculate-predicted-cells tm active-basal-segments active-apical-segments))
        (tm-active-basal-segments-set!    tm active-basal-segments)
        (tm-active-apical-segments-set!   tm active-apical-segments)
        (tm-matching-basal-segments-set!  tm matching-basal-segments)
        (tm-matching-apical-segments-set! tm matching-apical-segments)))))
                                                                                            ;
(define (activate-cells                  ;; TM {ColX} {CellX} {CellX} {CellX} {CellX} {CellX} {CellX} Boolean . {ColX} ->
          tm feedforward-input
          basal-reinforce-candidates apical-reinforce-candidates
          cc-basal-growth-candidates mc-basal-growth-candidates
          cc-apical-growth-candidates mc-apical-growth-candidates
          learning . bursting-columns)
  ;; update active-cells:      correct-predicted-cells + bursting columns
  ;;   winner-cells:           correct-predicted-cells + cells of learning-matching-basal-segments + new-basal-segment-cells
  ;;   predicted-active-cells: correct-predicted-cells
  ;; update synapses on active and matching segments, and create new segments
  (let* (                                ;; calculate active cells
      (correct-predicted-cells
        (cells-in-cols tm (tm-predicted-cells tm) feedforward-input))
      (bursting-columns
        (if (pair? bursting-columns)
          (car bursting-columns)
          (setdiff1d feedforward-input (cols-from-cells tm (tm-predicted-cells tm)))))
      (new-active-cells
        (append correct-predicted-cells (get-all-cells-in-columns tm bursting-columns))))
  (let*-values (                         ;; calculate learning
        [ ( learning-active-basal-segments
            learning-matching-basal-segments
            basal-segments-to-punish
            new-basal-segment-cells
            learning-cells)
              (calculate-basal-learning tm
                feedforward-input bursting-columns correct-predicted-cells) ]
        [ ( learning-active-apical-segments
            learning-matching-apical-segments
            apical-segments-to-punish
            new-apical-segment-cells)
              (calculate-apical-learning tm learning-cells feedforward-input) ] )
      (when learning
        (let ([basal-reinforce-candidates  (list->vector basal-reinforce-candidates)]
              [apical-reinforce-candidates (list->vector apical-reinforce-candidates)])
                                         ;; Learn on existing segments
          (learn tm learning-active-basal-segments basal-reinforce-candidates
              cc-basal-growth-candidates mc-basal-growth-candidates
              (tm-basal-pre-index tm) (tm-sample-size tm))
          (learn tm learning-matching-basal-segments basal-reinforce-candidates
              cc-basal-growth-candidates mc-basal-growth-candidates
              (tm-basal-pre-index tm) (tm-sample-size tm))
          (learn tm learning-active-apical-segments apical-reinforce-candidates
              cc-apical-growth-candidates mc-apical-growth-candidates
              (tm-apical-pre-index tm) (tm-apical-sample-size tm))
          (learn tm learning-matching-apical-segments apical-reinforce-candidates
              cc-apical-growth-candidates mc-apical-growth-candidates
              (tm-apical-pre-index tm) (tm-apical-sample-size tm))
                                         ;; Punish incorrect predictions
          (unless (zero? (tm-basal-predicted-segment-decrement tm))
            (adjust-active-synapses tm basal-segments-to-punish basal-reinforce-candidates
                                    (- (tm-basal-predicted-segment-decrement tm))))
          (unless (zero? (tm-apical-predicted-segment-decrement tm))
            (adjust-active-synapses tm apical-segments-to-punish apical-reinforce-candidates
                                    (- (tm-apical-predicted-segment-decrement tm)))))
                                         ;; Grow new segments
        (when (positive? (+ (length cc-basal-growth-candidates) (length mc-basal-growth-candidates)))
          (let ([new-segments
                  (create-segments tm (tm-basal-connections tm) new-basal-segment-cells)])
            (learn-on-new-segments tm new-segments
                cc-basal-growth-candidates mc-basal-growth-candidates
                (tm-basal-pre-index tm) (tm-sample-size tm))

            (for-each (lambda (cellx)
                (when (fxzero? (synapses-length (seg-synapses (car (vector-ref (tm-basal-connections tm) cellx)))))
                  (vector-set! (tm-basal-connections tm) cellx (cdr (vector-ref (tm-basal-connections tm) cellx)))))
              new-basal-segment-cells)))

        (when (positive? (+ (length cc-apical-growth-candidates) (length mc-apical-growth-candidates)))
          (let ([new-segments
                  (create-segments tm (tm-apical-connections tm) new-apical-segment-cells)])
            (learn-on-new-segments tm new-segments
                cc-apical-growth-candidates mc-apical-growth-candidates
                (tm-apical-pre-index tm) (tm-apical-sample-size tm)))

            (for-each (lambda (cellx)
                (when (fxzero? (synapses-length (seg-synapses (car (vector-ref (tm-apical-connections tm) cellx)))))
                  (vector-set! (tm-apical-connections tm) cellx (cdr (vector-ref (tm-apical-connections tm) cellx)))))
              new-apical-segment-cells)))

                                         ;; Save the results
      (tm-active-cells-set! tm (sort! fx<? new-active-cells))
      (tm-winner-cells-set! tm learning-cells)
      (tm-predicted-active-cells-set! tm correct-predicted-cells))))
                                                                                            ;
(define (calculate-basal-learning tm     ;; TM {ColX} {ColX} {CellX} -> {Seg} {Seg} {Seg} {CellX} {CellX} 
          feedforward-input bursting-columns correct-predicted-cells)
  ;; "Basic Temporal Memory" - see comments in .py
  (let* (
      (learning-active-basal-segments
        (filter-segments-by-cell (tm-active-basal-segments tm) correct-predicted-cells))
      (cells-for-matching-basal
        (map seg-cellx (tm-matching-basal-segments tm)))
      (basal-segments-to-punish          ;; here because unique! will mutate cells-for-matching-basal
        (exclude-segments tm (tm-matching-basal-segments tm) cells-for-matching-basal feedforward-input))
      (matching-cells (unique! fx=? cells-for-matching-basal))
      (matching-cells-in-bursting-columns
        (cells-in-cols tm matching-cells bursting-columns))
      (bursting-columns-with-no-match
        (setdiff1d bursting-columns (cols-from-cells tm matching-cells)))
      (learning-matching-basal-segments
        (choose-best-segment-per-column tm matching-cells-in-bursting-columns
          (tm-matching-basal-segments tm)))
      (new-basal-segment-cells
        (get-cells-with-fewest-segments tm bursting-columns-with-no-match))
      (learning-cells
        (unique! fx=? (sort! fx<?
          (append correct-predicted-cells
                  (map-segments-to-cells learning-matching-basal-segments)
                  new-basal-segment-cells
                  '())))))               ;; *copy* new-basal-segment-cells
    (values
      learning-active-basal-segments
      learning-matching-basal-segments
      basal-segments-to-punish
      new-basal-segment-cells
      learning-cells)))
                                                                                            ;
(define (calculate-apical-learning tm    ;; TM {CellX} {ColX} -> {Seg} {Seg} {Seg} {CellX}
          learning-cells feedforward-input)
  ;; calculate apical learning for each learning cell - see comments in .py
  (let* (
      (learning-active-apical-segments
        (filter-segments-by-cell (tm-active-apical-segments tm) learning-cells))
      (learning-cells-without-active-apical
        (setdiff1d learning-cells (unique! fx=? (map-segments-to-cells learning-active-apical-segments))))
      (cells-for-matching-apical
        (map-segments-to-cells (tm-matching-apical-segments tm)))
      (learning-cells-with-matching-apical
        (intersect1d learning-cells-without-active-apical cells-for-matching-apical))
      (learning-matching-apical-segments
        (choose-best-segment-per-cell learning-cells-with-matching-apical
                                      (tm-matching-apical-segments tm)))
      (new-apical-segment-cells
        (setdiff1d learning-cells-without-active-apical learning-cells-with-matching-apical))
      (apical-segments-to-punish
        (exclude-segments tm (tm-matching-apical-segments tm) cells-for-matching-apical feedforward-input)))
    (values
      learning-active-apical-segments
      learning-matching-apical-segments
      apical-segments-to-punish
      new-apical-segment-cells)))
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
(define (learn tm learning-segments      ;; TM  {Seg} (Vectorof CellX) {CellX} PreIndex ->
          active-input cc-growth-candidates mc-growth-candidates pre-index ss)
  ;; adjust synapse permanences, grow new synapses, and grow new segments
  (adjust-synapses tm learning-segments active-input)
  (let* ( 
      [max-new                           ;; Nat | {Nat}
        (if (fxnegative? ss)
          (fx+ (length cc-growth-candidates) (length mc-growth-candidates))
          (map (lambda (seg)
              (fx- ss (seg-potential-overlaps seg)))
            learning-segments))]
      [msps (tm-max-synapses-per-segment tm)]
      [max-new
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
                learning-segments max-new)))])
    (grow-synapses-to-sample tm learning-segments cc-growth-candidates mc-growth-candidates pre-index max-new)))
                                                                                            ;
(define (learn-on-new-segments tm        ;; TM Connections {CellX} {CellX} PreIndex ->
          new-segments cc-growth-candidates mc-growth-candidates pre-index ss)
  ;; connect new-segments to growth-candidates
  (let* ( 
      [num-new-synapses (fx+ (length cc-growth-candidates) (length mc-growth-candidates))]
      [num-new-synapses (if (fxnegative? ss)
                          num-new-synapses
                          (fxmin num-new-synapses ss))]
      [num-new-synapses (if (fxnegative? (tm-max-synapses-per-segment tm))
                          num-new-synapses
                          (fxmin num-new-synapses (tm-max-synapses-per-segment tm)))])
    (grow-synapses-to-sample tm new-segments cc-growth-candidates mc-growth-candidates pre-index num-new-synapses)))
                                                                                            ;
(define (choose-best-segment-per-cell    ;; {CellX} {Seg} Overlaps -> {Seg}
          cells all-matching-segments)
  ;; choose matching segment with max active potential synapses (on tie choose first)
  (let next-cell ([cells cells] [learning-segments '()])
    (if (null? cells) 
      (reverse! learning-segments)
      (let next-segment ( [segments all-matching-segments] 
                          [best-seg #f] 
                          [max-overlap (least-fixnum)])
        (cond
          [ (null? segments) (next-cell (cdr cells) (cons best-seg learning-segments)) ]
          [ (fx=? (car cells) (seg-cellx (car segments)))
              (let ([overlap (seg-potential-overlaps (car segments))])
                (if (fx>? overlap max-overlap)
                    (next-segment (cdr segments) (car segments) overlap)
                    (next-segment (cdr segments) best-seg max-overlap))) ]
          [ else (next-segment (cdr segments) best-seg max-overlap) ] )))))
                                                                                            ;
(define (choose-best-segment-per-column  ;; TM {CellX} {Seg} Overlaps -> {Seg}
          tm matching-cells all-matching-segments)
  ;; like choose-best-segment-per-cell but for cells in col
  (let* ( 
      [candidate-segments
        (filter-segments-by-cell all-matching-segments matching-cells)]
      [cell-scores
        (map (lambda (s) (seg-potential-overlaps s)) candidate-segments)]
      [columns-for-candidates
        (map (/cpc tm) (map-segments-to-cells candidate-segments))]
      [one-per-column-filter 
        (argmax-multi cell-scores columns-for-candidates)]
      [learning-segments
        (vector->list (vector-refs (list->vector candidate-segments)
                                   (list->vector one-per-column-filter)))])
    learning-segments))
                                                                                            ;
(define (get-cells-with-fewest-segments  ;; TM {ColX} -> {CellX}
          tm columnxs)
  ;; produce cell for each col with fewest total basal segments; break ties randomly
  (map (lambda (colx)
      (let loop ([cellx (fx* colx (tm-cells-per-column tm))]
                 [candidate-cells '()]
                 [fewest-segs (greatest-fixnum)])
        (if (fx<? cellx (fx* (fx1+ colx) (tm-cells-per-column tm)))
          (let ((n-segs (length (vector-ref (tm-basal-connections tm) cellx))))
            (cond
              [ (fx<? n-segs fewest-segs)
                  (loop (fx1+ cellx) (list cellx) n-segs) ]
              [ (fx=? n-segs fewest-segs)
                  (loop (fx1+ cellx) (cons cellx candidate-cells) n-segs) ]
              [ else (loop (fx1+ cellx) candidate-cells fewest-segs) ] ))
          (list-ref candidate-cells (random (length candidate-cells))))))
    columnxs))
    
;; === Supporting Functions ===
                                                                                            ;
(define (calculate-segment-activity tm   ;; TM {InpX} PreIndex {InpX} {CellX} -> {Seg} {Seg}
          active-input pre-index reduced-threshold-cells activation-threshold)
  ;; for each input index (may be external input or cell index) use pre-index (separate
  ;; indexes for basal/apical inputs) to find segments containing synapse for that input;
  ;; find synapse, build segment lists and update overlap counts according to thresholds
  ;; (combines attm.py calculateApical/BasalSegmentActivity and connections.py computeActivity)
  (define act-incr #x100)
  (define (act-count ov) (fxand #xFF00 ov))
  (define (pot-count ov) (fxand #x00FF ov))
  (let ([iteration (tm-iteration tm)]    ;; use iteration to clear previous overlap counts
        [actthresh (fx* act-incr activation-threshold)]
        [redthresh (if (null? reduced-threshold-cells)  -1
                       (fx* act-incr (tm-reduced-basal-threshold tm)))]
        [r-t-cells (list->vector reduced-threshold-cells)])
    (let inp-loop ([inpxs active-input] [psegs '()] [asegs '()])
      (if (pair? inpxs)
        (let ([inpx (car inpxs)])
          (let seg-loop ([segs (vector-ref pre-index inpx)] [psegs psegs] [asegs asegs])
            (if (pair? segs)
              (let* ( [segment (car segs)]
                      [synapse (synapses-search segment inpx)])
                (if synapse
                  (let ([new-overlap                                     ;; bump pot-count
                          (fx1+ (fxmax iteration (seg-overlap segment)))])
                    (seg-loop
                      (cdr segs)
                      (if (fx=? (pot-count new-overlap) (tm-min-threshold tm))
                        (cons segment psegs)
                        psegs)
                      (if (fx>=? (syn-perm synapse) (tm-connected-permanence tm))
                        (let ([new-overlap (fx+ new-overlap act-incr)])  ;; bump act-count
                          (seg-overlap-set! segment new-overlap)
                          (let ([overact (act-count new-overlap)])
                            (if (and (or (fx=? overact actthresh)
                                         (and (fx=? overact redthresh)
                                              (fxsearch r-t-cells (seg-cellx segment))))
                                     (not (memq segment asegs)))
                                (cons segment asegs)
                                asegs)))
                        (begin (seg-overlap-set! segment new-overlap)
                                asegs))))
                  (seg-loop (cdr segs) psegs asegs)))
              (inp-loop (cdr inpxs) psegs asegs))))
        (values
          asegs
          (unique! eq?
            (sort! (lambda (sega segb)
                (fx<=? (seg-cellx sega) (seg-cellx segb)))
              psegs)))))))
                                                                                            ;
(define (map-segments-to-cells segments) ;; {Seg} -> {CellX}
  ;; produce sorted cellxs from segments
  (sort! fx<? (map seg-cellx segments)))
                                                                                            ;
(define (filter-segments-by-cell         ;; {Seg} {CellX} -> {Seg}
          segments cellxs)
  ;; produce subset of segments that are on cellxs
  (filter (lambda (segment)
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
  (let per-group ([a a] [group-keys group-keys] [index 0] [result '()])
    (cond
      [ (null? group-keys) (reverse! result) ]
      [ else
        (let ([group (car group-keys)])
          (let next ([a a] [group-keys group-keys] [index index]
                     [max-in-group (least-fixnum)] [index-of-max #f])
            (cond
              [ (null? group-keys) (per-group a group-keys index (cons index-of-max result))]
              [ (fx=? (car group-keys) group)
                  (if (fx>? (car a) max-in-group)
                    (next (cdr a) (cdr group-keys) (fx1+ index) (car a) index)
                    (next (cdr a) (cdr group-keys) (fx1+ index) max-in-group index-of-max)) ]
              [ else (per-group a group-keys index (cons index-of-max result)) ]))) ])))
                                                                                            ;
(define (get-all-cells-in-columns tm     ;; TM {ColX} -> {CellX}
           colxs)
  ;; produce all cellxs in the colxs
  (let ([cpc (tm-cells-per-column tm)])
    (apply append! 
      (map (lambda (colx)
          (let ([first (fx* colx cpc)])
            (build-list cpc (lambda (i) (fx+ first i)))))
        colxs))))
                                                                                            ;
(define (prune-synapses tm segment omits);; Seg {Nat} ->
  ;; update segment, omitting synapses indexed by omits (which is sorted)
  (let ([syns (seg-synapses segment)])
    (let loop ([sx 0] [omits omits] [result (list)])
      (cond
        [ (fx=? sx (synapses-length syns))
            (seg-synapses<-list segment (reverse! result)) ]
        [ (if (null? omits)  #f  (fx=? sx (car omits)))
            (let ([prex (syn-prex (synapses-ref syns sx))])
              (remq! segment (vector-ref (tm-basal-pre-index tm) prex))
              (remq! segment (vector-ref (tm-apical-pre-index tm) prex)))
            (loop (fx1+ sx) (cdr omits) result) ]
        [ else (loop (fx1+ sx) omits (cons (synapses-ref syns sx) result)) ]))))
                                                                                            ;
(define (adapt-segment tm segment        ;; TM Seg (Vectorof CellX) Perm Perm ->
          active-input permanence-delta permanence-decrement)
  ;; Updates synapses on segment. Strengthens active synapses; weakens inactive synapses.
  ;; Remove synapse on zero permanence, destroy segment if no synapses left.
  (let ([synapses (seg-synapses segment)]
        [trim-permanence (tm-perm-trim-threshold tm)])
    (let build-s-t-d ([sx (fx1- (synapses-length synapses))]
                      [synapses-to-destroy '()])
      (if (fxnegative? sx)
          (cond [ (null? synapses-to-destroy) ]
                [ (fx=? (length synapses-to-destroy) (synapses-length synapses))
                    (for-each (lambda (s)
                        (remq! segment (vector-ref (tm-basal-pre-index tm) (syn-prex s)))
                        (remq! segment (vector-ref (tm-apical-pre-index tm) (syn-prex s))))
                      (seg-synapses->list segment))
                    (seg-synapses<-list segment '()) ]
                [ else (prune-synapses tm segment synapses-to-destroy) ])
          (let* ( [synapse (synapses-ref synapses sx)]
                  [prex    (syn-prex synapse)]
                  [perm    (if (fxsearch active-input prex)
                             (clip-max (clip-min (fx+ (syn-perm synapse) permanence-delta)))
                             (clip-min (fx- (syn-perm synapse) permanence-decrement)))])
            (if (fx<=? perm trim-permanence)
                ;; build synapses-to-destroy indices as sorted list
                (build-s-t-d (fx1- sx) (cons sx synapses-to-destroy))
                (begin
                  (synapses-set! synapses sx (make-syn prex perm))
                  (build-s-t-d (fx1- sx) synapses-to-destroy))))))))
                                                                                            ;
(define (adjust-synapses tm segments     ;; TM {Seg} (Vectorof CellX) ->
          active-input)
  ;; For each specified segment, update the permanence of each synapse
  ;; according to whether the synapse would be active given the specified
  ;; active inputs.
  (let ([inc (tm-permanence-increment tm)]
        [dec (tm-permanence-decrement tm)])
    (for-each (lambda (segment)
        (adapt-segment tm segment active-input inc dec))
      segments)))
                                                                                            ;
(define (adjust-active-synapses tm       ;; TM {Seg} (Vectorof CellX) Perm ->
          segments active-input permanence-delta)
  ;; For each specified segment, add a delta to the permanences of the
  ;; synapses that would be active given the specified active inputs.
  (for-each (lambda (segment)
      (adapt-segment tm segment active-input permanence-delta 0))
    segments))
                                                                                            ;
(define (grow-synapses-to-sample tm      ;; TM {Seg} {CellX} {CellX} PreIndex Nat|{Nat} ->
          segments cc-inputs mc-inputs pre-index sample-size)
  ;; For each segment, grow synapses to a random subset of the
  ;; inputs that aren't already connected to the segment.
  (if (fixnum? sample-size)
    (for-each (lambda (segment)
        (grow-synapses tm segment cc-inputs mc-inputs pre-index sample-size))
      segments)
    (for-each (lambda (segment sample-size)
        (grow-synapses tm segment cc-inputs mc-inputs pre-index sample-size))
      segments sample-size)))
                                                                                            ;
(define (add-to-pre-index                ;; CellX Seg PreIndex ->
          prex seg pre-index)
  ;; add to the list of segments with synapses from the pre-synaptic cell
  (let ((pre-list (vector-ref pre-index prex)))
    (unless (memq seg pre-list)
      (vector-set! pre-index prex (cons seg pre-list)))))
                                                                                            ;
(define (grow-synapses                   ;; TM Seg {CellX} {CellX} PreIndex Nat ->
          tm segment cc-inputs mc-inputs pre-index n-desired-new-synapses)
  ;; grow synapses to all inputs that aren't already connected to the segment
  (let build-candidates ([cis '()] [ccis cc-inputs] [mcis mc-inputs] [nc 0])  ;; candidate inputs
    (cond 
      [ (and (null? ccis) (null? mcis))
        (if (fx<=? nc n-desired-new-synapses)
          (let ([init-perm (tm-initial-permanence tm)])
            (seg-synapses<-list segment 
              (merge! fx<?
                (sort! fx<?
                  (let map-car! ([cic cis])   ;; candidate inputs cursor
                    (cond [(null? cic) cis]
                      [else
                        (add-to-pre-index (car cic) segment pre-index)
                        (set-car! cic (make-syn (car cic) init-perm))
                        (map-car! (cdr cic)) ])))
                (seg-synapses->list segment))))
          (when (fxpositive? n-desired-new-synapses)
            (let ([new-prexs (vector-sample (list->vector cis) n-desired-new-synapses)]
                  [init-perm (tm-initial-permanence tm)])
              (vector-sort! fx<? new-prexs)
              (seg-synapses<-list segment 
                (merge! fx<? 
                  (build-list n-desired-new-synapses (lambda (newx) 
                      (let ((new-prex (vector-ref new-prexs newx)))
                        (add-to-pre-index new-prex segment pre-index)
                        (make-syn new-prex init-perm))))
                  (seg-synapses->list segment)))))) ]
      [ (pair? ccis) 
          (if (synapses-search segment (car ccis))
              (build-candidates cis (cdr ccis) mcis nc)
              (build-candidates (cons (car ccis) cis) (cdr ccis) mcis (fx1+ nc))) ]
      [ else
          (if ((tm-same-col? tm) (fxdiv (seg-cellx segment) (tm-cells-per-column tm)) (car mcis))
            (if (synapses-search segment (car mcis))
                (build-candidates cis ccis (cdr mcis) nc)
                (build-candidates (cons (car mcis) cis) ccis (cdr mcis) (fx1+ nc)))
            (build-candidates cis ccis (cdr mcis) nc)) ])))
                                                                                            ;
(define (create-segment tm               ;; TM Connections CellX -> Seg
          connections cellx)
  ;; produce a new segment on the cell
  (let ([segments (vector-ref connections cellx)]
        [segment (make-seg cellx)])
    (vector-set! connections cellx (cons segment segments))
    segment))
                                                                                            ;
(define (create-segments tm              ;; TM Connections {CellX} -> {Seg}
          connections cellxs)
  (map (lambda (cellx)
      (create-segment tm connections cellx))
    cellxs))
                                                                                            ;
(define (number-of-cells tm)             ;; TM -> Nat
  (fx* (tm-column-count tm) (tm-cells-per-column tm)))
                                                                                            ;
(define (number-of-segments connections) ;; (CellVecOf {Seg}) -> Nat
  (vector-fold-left (lambda (acc segments)
      (+ acc (length segments)))
    0
    connections))
                                                                                            ;
(define (number-of-synapses connections) ;; (CellVecOf {Seg}) -> Nat
  (vector-fold-left (lambda (acc segments)
      (+ acc
        (fold-left (lambda (acc segment)
            (+ acc (synapses-length (seg-synapses segment))))
          0
          segments)))
    0
    connections))
                                                                                            ;
(define (number-of-basal-segments tm)    ;; TM -> Nat
  (number-of-segments (tm-basal-connections tm)))
                                                                                            ;
(define (number-of-apical-segments tm)   ;; TM -> Nat
  (number-of-segments (tm-apical-connections tm)))
                                                                                            ;
(define (number-of-basal-synapses tm)    ;; TM -> Nat
  (number-of-synapses (tm-basal-connections tm)))
                                                                                            ;
(define (number-of-apical-synapses tm)   ;; TM -> Nat
  (number-of-synapses (tm-apical-connections tm)))
                                                                                            ;
(define (average-preindex-length tm      ;; TM PreIndex -> Number
          pre-index)
  (let ([counts
    (vector-fold-left (lambda (acc segs)
        (cons (+ (car acc) (length segs))
              (+ (cdr acc) (if (null? segs) 0 1))))
      '(0 . 0)
      (pre-index tm))])
    (if (zero? (cdr counts))  0
      (/ (inexact (quotient (* 10 (car counts)) (cdr counts))) 10))))

)