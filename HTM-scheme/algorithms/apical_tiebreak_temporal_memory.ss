;; HTM-scheme Apical Tiebreak Temporal Memory algorithm (C) 2019-2021 Roger Turner.
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
  activate-cells has an optional bursting-columns input (see layer4.ss for usage)
  a connect? hook in grow-synapses can be used to filter growth-candidates
  winner-cells in .py renamed learning-cells

*See htm_concept.ss for type and data structure descriptions and code conventions*
  
Description from apical_tiebreak_temporal_memory.py:
A generalized Temporal Memory with apical dendrites that add a "tiebreak".
Basal connections are used to implement traditional Temporal Memory.
The apical connections are used for further disambiguation. If multiple cells
in a minicolumn have active basal segments, each of those cells is predicted,
unless one of them also has an active apical segment, in which case only the
cells with active basal and apical segments are predicted.

"Remember that all models are wrong; the practical question is how wrong
do they have to be to not be useful." [George Box]
  
  |#

  #!chezscheme

(library (HTM-scheme HTM-scheme algorithms apical_tiebreak_temporal_memory)
                                                                                            ;
(export
make-tm
reset
depolarize-cells
activate-cells
  (rename
    (tm-active-cells           get-active-cells)
    (tm-learning-cells         get-learning-cells)
    (tm-predicted-active-cells get-predicted-active-cells)
    (tm-predicted-cells        get-predicted-cells))
  tm                                     ;; for pair/sequence memory
  tm-axon-radius2
  tm-connected-permanence
  number-of-cells
  number-of-connected-cells
  number-of-basal-segments
  number-of-basal-synapses
  number-of-apical-segments
  number-of-apical-synapses
  connection-lengths
  get-axon-tree
  cellx->colx
  map-segments-to-cells
  projection-of-source
  cells-in-cols
  cols-from-cells
  )
                                                                                            ;
(import
  (rename (chezscheme) (sort! sort-by!) (reset $reset))
  (parameters)
  (HTM-scheme HTM-scheme algorithms htm_prelude)
  (HTM-scheme HTM-scheme algorithms htm_concept)
  (HTM-scheme HTM-scheme math       coordinates))
                                                                                            ;
  (implicit-exports #f)

;; === Types (see htm_concept.ss for common HTM types) ===
                                                                                            ;
;; Dendrite     = Symbol: 'proximal | 'basal | 'apical
;; TM           = Record: tm parameters
;; PBA          = Record: proximal/basal/apical parameters

;; === Layer data structures ===
                                                                                            ;
(define-record-type tm (fields           ;; TM [Temporal Memory parameters: cf __init__ in attm.py]
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
    max-segments-per-cell                ;; Fixnum [-1 => unlimited]
    seed                                 ;; Integer
    (mutable active-cells)               ;; (Listof CellX)
    (mutable learning-cells)             ;; (Listof CellX)
    (mutable predicted-cells)            ;; (Listof CellX)
    (mutable predicted-active-cells)     ;; (Listof CellX)
    proximal-ds                          ;; PBA [proximal dendrite parameters]
    basal-ds                             ;; PBA [basal dendrite parameters]
    apical-ds                            ;; PBA [apical dendrite parameters]
    use-apical-tiebreak                  ;; Boolean
    use-apical-modulation-basal-threshold;; Boolean
    ;; following parameters are not in apical_tiebreak_temporal_memory.py
    (mutable iteration)                  ;; Nat [incrementing by #x10000]
    cortical-column                      ;; Nat
    layer                                ;; Nat
    axon-radius2                         ;; Nat
    perm-trim-threshold                  ;; Perm
    use-bursting-columns-input           ;; Boolean
    basal-connect?                       ;; (CellX ColX -> Boolean) | #f
    apical-connect?                      ;; (CellX ColX -> Boolean) | #f
    number-of-basal-segments)            ;; (FXVector CellX->Nat)   | #f
  (opaque #t) (nongenerative tm)
(protocol #;(make-tm kwargs)             ;; {KWarg} -> TM
  ;; use kwargs to override defaults, construct initial values, create record
  (lambda (new)
    (lambda (kwargs)
      (let* (
          [tm       (apply new (key-word-args kwargs attm-defaults))]
          [n-cells  (fx* (tm-column-count tm) (tm-cells-per-column tm))]
          [kwargs   (if (tm-number-of-basal-segments tm)
                      (cons `[number-of-basal-segments . ,(make-fxvector n-cells 0)] kwargs)
                      kwargs)]
          [kwargs   (cons `[basal-ds . ,(make-pba
                              (tm-basal-input-size tm)
                              (tm-basal-predicted-segment-decrement tm)
                              (make-at n-cells)
                              (tm-basal-connect? tm)
                              (make-eqv-hashtable n-cells)
                              (list)
                              (list))] kwargs)]
          [kwargs   (cons `[apical-ds . ,(make-pba
                              (tm-apical-input-size tm)
                              (tm-apical-predicted-segment-decrement tm)
                              (make-at n-cells)
                              (tm-apical-connect? tm)
                              (make-eqv-hashtable n-cells)
                              (list)
                              (list))] kwargs)])
        (apply new (key-word-args kwargs attm-defaults)))))))
                                                                                            ;
(define attm-defaults `(                 ;; (Listof KWarg) [cf __init__ in attm.py]
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
    [max-segments-per-cell                 . 255]
    [seed                                  . 42]
    [active-cells                          . ()]
    [learning-cells                        . ()]
    [predicted-cells                       . ()]
    [predicted-active-cells                . ()]
    [proximal-ds                           . #f]
    [basal-ds                              . #f]
    [apical-ds                             . #f]
    [use-apical-tiebreak                   . #t]
    [use-apical-modulation-basal-threshold . #t]
    [iteration                             . 0]
    [cortical-column                       . 0]
    [layer                                 . 0]
    [axon-radius2                          . 1]
    [perm-trim-threshold                   . ,(perm 0.0)]  ;; default value => attm.py behaviour
    [use-bursting-columns-input            . #f]           ;; ''
    [basal-connect?                        . #f]           ;; ''
    [apical-connect?                       . #f]           ;; ''
    [number-of-basal-segments              . #vfx()] ))    ;; ''
                                                                                            ;
(define-record-type pba (fields          ;; PBA [proximal/basal/apical parameters]
    ;; initialised in make-tm; pos-cache caches column lists for projection-of-source
    input-size                           ;; Nat
    predicted-segment-decrement          ;; Perm
    axon-tree                            ;; AxonTree
    connect?                             ;; (Source ColX -> Boolean) | #f
    pos-cache                            ;; (Hashtable Source->[Nat . {ColX}])
    (mutable active-segments)            ;; {Segment}
    (mutable matching-segments))         ;; {Segment}
  (sealed #t) (opaque #t) (nongenerative pba))
                                                                                            ;
(define (ba dend)                        ;; Dendrite -> PBA
  ;; produce basal-ds or apical-ds 
  (case dend
    [(basal) tm-basal-ds]
    [else    tm-apical-ds]))

;; === Apical Tiebreak Temporal Memory Algorithm ===
                                                                                            ;
(define (reset tm)                       ;; TM ->
  ;; Clear all cell and segment activity (overlap counts are stored in segments)
  (tm-active-cells-set!              tm (list))
  (tm-learning-cells-set!            tm (list))  ;; htmresearch "winner" renamed learning
  (tm-predicted-cells-set!           tm (list))
  (tm-predicted-active-cells-set!    tm (list))
  (pba-active-segments-set!   (tm-basal-ds  tm) (list))
  (pba-active-segments-set!   (tm-apical-ds tm) (list))
  (pba-matching-segments-set! (tm-basal-ds  tm) (list))
  (pba-matching-segments-set! (tm-apical-ds tm) (list)))
                                                                                            ;
(define (depolarize-cells tm             ;; TM {CellX} {CellX} Boolean ->
          basal-input                    ;; active cells pre-synaptic to basal segments
          apical-input                   ;; active cells pre-synaptic to apical segments
          learn)                         ;; #f (inference mode) to use reduced-basal-threshold
  ;; Calculate predictions. Save with active/matching segments (overlaps saved in segments)
  (tm-iteration-set! tm (fx+ iteration-incr (tm-iteration tm)))
                                         ;; (see calculate-segment-activity for iteration-incr use)
  (when (fxpositive? (tm-cells-per-column tm))
    (let*-values (
        [(active-apical-segments matching-apical-segments)
           (calculate-segment-activity
             tm apical-input '() (pba-axon-tree (tm-apical-ds tm)))]
        [(reduced-basal-threshold-cells)
           (if (or learn (not (tm-use-apical-modulation-basal-threshold tm)))  '()
               (map-segments-to-cells active-apical-segments))]
        [(active-basal-segments matching-basal-segments)
           (calculate-segment-activity
             tm basal-input reduced-basal-threshold-cells (pba-axon-tree (tm-basal-ds tm)))])
      (tm-predicted-cells-set! tm
        (calculate-predicted-cells tm active-basal-segments active-apical-segments))
      (pba-active-segments-set!   (tm-basal-ds  tm) active-basal-segments)
      (pba-active-segments-set!   (tm-apical-ds tm) active-apical-segments)
      (pba-matching-segments-set! (tm-basal-ds  tm) matching-basal-segments)
      (pba-matching-segments-set! (tm-apical-ds tm) matching-apical-segments))))
                                                                                            ;
(define (activate-cells tm               ;; TM {ColX} {CellX} {CellX} {CellX} {CellX} Boolean . {ColX} ->
          active-columns                 ;; minicolumns with proximal input
          basal-reinforce-candidates     ;; pre-synaptic cells whose basal synapses should be incremented
          apical-reinforce-candidates    ;; pre-synaptic cells whose apical synapses should be incremented
          basal-growth-candidates        ;; potential pre-synaptic cells for basal segments
          apical-growth-candidates       ;; potential pre-synaptic cells for apical segments
          learn . bursting-columns)      ;; update segments/synapses? optional bursting columns input
  ;; Activate cells in the specified columns, using predictions from depolarize-cells, then learn.
  ;;   active-cells:           correct-predicted-cells + all cells in bursting columns
  ;;   learning-cells:         correct-predicted-cells + cells of learning-matching-basal-segments + new-basal-segment-cells
  ;;   predicted-active-cells: correct-predicted-cells
  ;; update synapses on active and matching segments, and create new segments
#;> (define (learn-pba pba               ;; PBA {Seg} {Seg} {Seg} {CellX} {Source} {Source} ->
              learning-active-segments learning-matching-segments
              segments-to-punish new-segment-cells
              reinforce-candidates growth-candidates)
      ;; Adapt learning segments, punish incorrect predictions, grow new segments
#;>   (define (create-segment cellx)     ;; CellX -> Segment
        ;; Produce new segment for the cell; count basal segments for use by get-cells-with-fewest-segments
        (let ([nobs (tm-number-of-basal-segments tm)])
          (when (and nobs (eq? pba (tm-basal-ds tm)))
            (fxvector-set! nobs cellx (fx1+ (fxvector-ref nobs cellx)))))
        (at-make-seg (pba-axon-tree pba) (tm-cortical-column tm) cellx))
      ;; (learn-pba)
      (let ([learning-segments    (append! learning-active-segments learning-matching-segments)]
            [reinforce-candidates (list->fxvector reinforce-candidates)]
            [increment            (tm-permanence-increment tm)]
            [decrement            (tm-permanence-decrement tm)])
        (for-each (lambda (segment)        ;; Learn on existing segments
            (adapt-segment tm segment reinforce-candidates increment decrement pba))
          learning-segments)
        (let ([decrement (fx- (pba-predicted-segment-decrement pba))])
          (unless (fxzero? decrement)
            (for-each (lambda (segment)    ;; Punish incorrect predictions
                (adapt-segment tm segment reinforce-candidates decrement 0 pba))
              segments-to-punish)))
        (grow-synapses-to-sample tm learning-segments growth-candidates pba)
        (when (pair? growth-candidates)    ;; Grow new segments
          (let ([new-segments (map create-segment new-segment-cells)])
            (grow-synapses-to-sample tm new-segments growth-candidates pba)))))
  ;; (activate-cells)                    ;; Calculate active cells
  (when (fxpositive? (tm-cells-per-column tm))
    (let* (
        [correct-predicted-cells
          (cells-in-cols tm (tm-predicted-cells tm) active-columns)]
        [bursting-columns
          (if (and (tm-use-bursting-columns-input tm) (pair? bursting-columns))
            (car bursting-columns)
            (setdiff1d active-columns (cols-from-cells tm (tm-predicted-cells tm))))]
        [new-active-cells
          (append correct-predicted-cells (all-cells-in-columns tm bursting-columns))])
      (if learn
        (let*-values (                   ;; Calculate learning
            [ ( learning-active-basal-segments
                learning-matching-basal-segments
                basal-segments-to-punish
                new-basal-segment-cells
                learning-cells)
                  (calculate-basal-learning tm
                    active-columns bursting-columns correct-predicted-cells) ]
            [ ( learning-active-apical-segments
                learning-matching-apical-segments
                apical-segments-to-punish
                new-apical-segment-cells)
                  (calculate-apical-learning tm learning-cells active-columns) ] )
          (learn-pba (tm-basal-ds tm) learning-active-basal-segments learning-matching-basal-segments
            basal-segments-to-punish new-basal-segment-cells
            basal-reinforce-candidates basal-growth-candidates)
          (learn-pba (tm-apical-ds tm) learning-active-apical-segments learning-matching-apical-segments
            apical-segments-to-punish new-apical-segment-cells
            apical-reinforce-candidates apical-growth-candidates)
                                         ;; Save the results
          (tm-learning-cells-set! tm learning-cells))
        ;; (not learn)
        (tm-learning-cells-set! tm '()))
      (tm-active-cells-set!           tm (sort-unique! new-active-cells))
      (tm-predicted-active-cells-set! tm correct-predicted-cells))))
                                                                                            ;
(define (calculate-basal-learning tm     ;; TM {ColX} {ColX} {CellX} -> {Seg} {Seg} {Seg} {CellX} {CellX} 
          active-columns bursting-columns correct-predicted-cells)
  ;; Basic Temporal Memory learning: see comments in .py
  (let* (
      [learning-active-basal-segments
        (filter-segments-by-cell (pba-active-segments (tm-basal-ds tm)) correct-predicted-cells)]
      [matching-basal-segments
        (pba-matching-segments (tm-basal-ds tm))]
      [cells-for-matching-basal
        (map-sorted-segments-to-cells matching-basal-segments)]
      [matching-cells-in-bursting-columns
        (cells-in-cols tm cells-for-matching-basal bursting-columns)]
      [bursting-columns-with-no-match
        (setdiff1d bursting-columns (cols-from-cells tm cells-for-matching-basal))]
      [learning-matching-basal-segments
        (choose-best-segment-per-column tm matching-basal-segments matching-cells-in-bursting-columns)]
      [new-basal-segment-cells
        (if (tm-number-of-basal-segments tm)
          (get-cells-with-fewest-segments tm bursting-columns-with-no-match)
          (random-cell-in-columns tm bursting-columns-with-no-match))]
      [learning-cells (sort-unique!
          (append correct-predicted-cells new-basal-segment-cells
                  (map seg-cellx learning-matching-basal-segments)))]
      [basal-segments-to-punish
        (exclude-segments tm matching-basal-segments active-columns)])
    (values
      learning-active-basal-segments
      learning-matching-basal-segments
      basal-segments-to-punish
      new-basal-segment-cells
      learning-cells)))
                                                                                            ;
(define (calculate-apical-learning tm    ;; TM {CellX} {ColX} -> {Seg} {Seg} {Seg} {CellX}
          learning-cells active-columns)
  ;; Calculate apical learning for each learning cell - see comments in .py
  (let* (
      [learning-active-apical-segments
        (filter-segments-by-cell (pba-active-segments (tm-apical-ds tm)) learning-cells)]
      [learning-cells-without-active-apical
        (setdiff1d learning-cells (map-segments-to-cells learning-active-apical-segments))]
      [matching-apical-segments
        (pba-matching-segments (tm-apical-ds tm))]
      [cells-for-matching-apical
        (map-sorted-segments-to-cells matching-apical-segments)]
      [learning-cells-with-matching-apical
        (intersect1d learning-cells-without-active-apical cells-for-matching-apical)]
      [learning-matching-apical-segments
        (choose-best-segment-per-cell tm learning-cells-with-matching-apical matching-apical-segments)]
      [new-apical-segment-cells
        (setdiff1d learning-cells-without-active-apical learning-cells-with-matching-apical)]
      [apical-segments-to-punish
        (exclude-segments tm matching-apical-segments active-columns)])
    (values
      learning-active-apical-segments
      learning-matching-apical-segments
      apical-segments-to-punish
      new-apical-segment-cells)))
                                                                                            ;
#;(calculate-apical-segment-activity)    ;; see calculate-segment-activity
#;(calculate-basal-segment-activity)     ;; see calculate-segment-activity
                                                                                            ;
(define (calculate-predicted-cells tm    ;; TM {Seg} {Seg} -> {CellX}
          active-basal-segments active-apical-segments)
  ;; Calculate predicted cells for given active segments
  (let ([cells-for-basal-segments (map-segments-to-cells active-basal-segments)])
    (if (not (tm-use-apical-tiebreak tm))
      cells-for-basal-segments
      (let* (
          [cells-for-apical-segments (map-segments-to-cells active-apical-segments)]
          [fully-depolarized-cells   (intersect1d cells-for-basal-segments
                                                  cells-for-apical-segments)]
          [partly-depolarized-cells  (setdiff1d cells-for-basal-segments
                                                fully-depolarized-cells)]
          [cellx->colx               (lambda (cellx) (cellx->colx tm cellx))]
          [inhibited-mask            (in1d (map cellx->colx partly-depolarized-cells)
                                           (map cellx->colx fully-depolarized-cells))])
        (sort-unique!
          (append! fully-depolarized-cells
                   (exclude-by-mask partly-depolarized-cells inhibited-mask)))))))
                                                                                            ;
#;(learn, learn-on-new-segments)         ;; see learn-pba in activate-cells
                                                                                            ;
(define (choose-best-segment-per-cell tm ;; TM {CellX} {Seg} -> {Seg}
          cells all-matching-segments)
  ;; Choose matching segment with max active potential synapses; previous code (cf attm.py):
  ;; (let* ([candidate-segments  (filter-segments-by-cell all-matching-segments cells)]
  ;;        [one-per-cell-filter (argmax-multi (map (lambda (s) (seg-potential-overlaps tm s)) candidate-segments)
  ;;                                           (map seg-cellx candidate-segments))])
  ;;        [learning-segments   (list-refs! candidate-segments one-per-cell-filter)])
  (filter-max-by-group tm all-matching-segments cells seg-cellx))
                                                                                            ;
(define (choose-best-segment-per-column  ;; TM {Seg} {CellX} -> {Seg}
          tm all-matching-segments matching-cells)
  ;; Like choose-best-segment-per-cell but for cells in col; previous code (cf attm.py):
  ;; (let* ([candidate-segments     (filter-segments-by-cell all-matching-segments matching-cells)]
  ;;        [cell-scores            (map (lambda (s) (seg-potential-overlaps tm s)) candidate-segments)]
  ;;        [columns-for-candidates (map (lambda (s) (cellx->colx tm (seg-cellx s))) candidate-segments)]
  ;;        [one-per-column-filter  (argmax-multi cell-scores columns-for-candidates)])
  ;;      (list-refs! candidate-segments one-per-column-filter))
  (filter-max-by-group tm all-matching-segments matching-cells
    (lambda (seg) (cellx->colx tm (seg-cellx seg)))))
                                                                                            ;
(define (get-cells-with-fewest-segments  ;; TM {ColX} -> {CellX}
          tm columnxs)
  ;; Produce cell for each col with fewest basal segments; break ties randomly
  (map (lambda (colx)
      (let* ( [start (fx* colx (tm-cells-per-column tm))]
              [limit (fx+ start (tm-cells-per-column tm))])
        (let loop ([cellx start] [cellxs (list)] [fewest (greatest-fixnum)])
          (if (fx<? cellx limit)
            (let ([n-segs (fxvector-ref (tm-number-of-basal-segments tm) cellx)]
                  [nextx  (fx1+ cellx)])
              (cond
                [ (fx<? n-segs fewest) (loop nextx (list cellx)        n-segs) ]
                [ (fx=? n-segs fewest) (loop nextx (cons cellx cellxs) n-segs) ]
                [ else                 (loop nextx cellxs              fewest) ] ))
            (if (null? (cdr cellxs))  (car cellxs)
                (list-ref cellxs (random (length cellxs))))))))
    columnxs))
                                                                                            ;
(define (random-cell-in-columns tm colxs);; TM {ColX} -> {CellX}
  ;; Produce random cell in each specified column (alternative to get-cells-with-fewest-segments)
  (let ([cpc (tm-cells-per-column tm)])
    (map (lambda (colx)
        (let ([first (fx* colx cpc)])
          (fx+ first (random cpc))))
      colxs)))

;; === Connections ===
                                                                                            ;
(define (adapt-segment tm segment        ;; TM Seg (FXVectorOf Source) Perm Perm PBA ->
          reinforce-candidates permanence-delta permanence-decrement pba)
  ;; Update segment's synapses: strengthen active, weaken inactive synapses;
  ;; remove synapse on low permanence, remove segment from indexes if few synapses left
  (define (remove-reference synapse)     ;; Synapse -> #f
    ;; Remove reference to segment from axon tree and clear pos-cache for this source
    (let* ( [axon-tree (pba-axon-tree pba)]
            [source    (syn-source synapse)]
            [segxv     (at-segxv axon-tree source)])
      (when segxv
        (segxv-remove! (lambda (segx)
            (eq? segment (at-seg-ref axon-tree segx)))
          segxv))
      (hashtable-delete! (pba-pos-cache pba) source))
    #f )
  ;; (adapt-segment)
  (let* (
      [trim-permanence (tm-perm-trim-threshold tm)]
      [synapses
        (synapses-update! (lambda (synapse)
            (let* ( [source (syn-source synapse)]
                    [perm   (if (fxsearch reinforce-candidates source)
                              (clip-max (fx+ (syn-perm synapse) permanence-delta))
                              (fx- (syn-perm synapse) permanence-decrement))])
              (if (fx<=? perm trim-permanence)
                (remove-reference synapse)   ;; synapses-update! will remove synapse
                (make-syn source perm))))
          (seg-synapses segment))])
    (when (fx>? (tm-min-threshold tm) (synapses-length synapses))
      (synapses-for-each remove-reference synapses)
      (seg-synapses-set! segment (make-synapses 0))
      (at-seg-set! (pba-axon-tree pba) (seg-segx segment) #f)
      (let ([nobs  (tm-number-of-basal-segments tm)]
            [cellx (seg-cellx segment)])
        (when (and nobs (eq? pba (tm-basal-ds tm)))
          (fxvector-set! nobs cellx (fx1- (fxvector-ref nobs cellx))))))))
                                                                                            ;
(define (grow-synapses-to-sample tm      ;; TM {Segment} {Source} PBA ->
          segments growth-candidates pba)
  ;; Add synapses (not already connected) for selection of growth-candidates to segments
  (define (partition-gc colx)            ;; ColX -> {Source} {Source}
    ;; divide growth-candidates into preferred (already connect near segment's col) and others
    (if no-connectivity?  (values '() growth-candidates)
        (let ([ccx (tm-cortical-column tm)])
          (partition (lambda (gc)
              (exists (lambda (projected-to-colx)
                  (fx>=? 1 (within-cc-distance2 colx projected-to-colx)))
                (let ([pos (projection-of-source tm gc pba)])
                  (if (fx=? (source-ccx gc) ccx)
                    (cons (cellx->colx tm (source-cellx gc)) pos)
                    pos))))
            growth-candidates))))
#;> (define (grow-synapses segment colx preferred others)
      ;; Add synapses that could connect
      (define (connectable colx candidates)
        ;; exclude candidates by connect rules, or already connected
        (let ([connect? (pba-connect? pba)]
              [synapses (seg-synapses segment)])
          (if connect?
            (filter (lambda (gc)
                (and (connect? gc colx) (not (synapses-search gc synapses))))
              candidates)
            (remp (lambda (gc) (synapses-search gc synapses)) candidates))))
      (let* (
          [preferred    (connectable -1 preferred)]
          [n-preferred  (length preferred)]
          [n-desired-new-synapses
            (if (fxnegative? (tm-sample-size tm))
              (greatest-fixnum)          ;; (potential-overlaps is zero for new segments)
              (fxmax 0 (fx- (tm-sample-size tm) (seg-potential-overlaps tm segment))))]
          [n-desired-new-synapses
            (if (fxnegative? (tm-max-synapses-per-segment tm))  n-desired-new-synapses
                (fxmin n-desired-new-synapses
                  (fx- (tm-max-synapses-per-segment tm) (synapses-length (seg-synapses segment)))))]
          [preferred    (if (fx<=? n-preferred n-desired-new-synapses)  preferred
                            (u32-sample preferred n-desired-new-synapses))]
          [n-others     (fx- n-desired-new-synapses n-preferred)]
          [others       (if (fxnonpositive? n-others)  '()
                            (u32-sample (connectable colx others) n-others))]
          [growth-candidates   (append preferred others)]
          [n-growth-candidates (length growth-candidates)])
        (when (fxpositive? n-growth-candidates)
          (let ([axon-tree (pba-axon-tree pba)])
            (for-each (lambda (gc)
                (at-update! axon-tree gc segment))
              growth-candidates))
          (seg-synapses-set! segment
            (synapses-add!
              (sort-unique-n!
                growth-candidates
                n-growth-candidates)
              n-growth-candidates
              (tm-initial-permanence tm)
              (seg-synapses segment))))))
  ;; (grow-synapses-to-sample): memoize preferred & others which depend only on segment's col
  (let next-seg ([segs segments] [colx-memo (list)])
    (when (pair? segs)
      (let* ( [segment (car segs)]
              [colx    (cellx->colx tm (seg-cellx segment))]
              [memo    (assv colx colx-memo)])
        (if memo
          (let ([memo (cdr memo)])
            (grow-synapses segment colx (car memo) (cdr memo))
            (next-seg (cdr segs) colx-memo))
          (let-values ([(preferred others) (partition-gc colx)])
            (grow-synapses segment colx preferred others)
            (next-seg (cdr segs)
              (cons (cons colx (cons preferred others)) colx-memo))))))))
                                                                                            ;
  ;; Pack overlap counts and iteration "time-stamp" into one segment field
  (define pot-overlap-mask #x000FF)
  (define act-overlap-mask #x0FF00)
  (define act-overlap-incr #x00100)
  (define iteration-incr   #x10000)
                                                                                            ;
(define (seg-potential-overlaps tm seg)  ;; TM Seg -> Nat
  ;; Produce potential overlaps for segment (ignoring value from a previous iteration)
  (let ([overlap (seg-overlap seg)])
    (if (fx<? overlap (tm-iteration tm))  0
        (fxand pot-overlap-mask overlap))))
                                                                                            ;
(define (calculate-segment-activity tm   ;; TM {Source} {CellX} AxonTree -> {Seg} {Seg}
          active-cells reduced-threshold-cells at)
  ;; For each source in active-cells use basal/apical axon tree to find segments with
  ;; synapses for that input, build segment lists and update overlap counts per thresholds
  ;; (combines attm.py calculateApical/BasalSegmentActivity and connections.py computeActivity)
  ;(define (csa-thunk)                    ;; (for cost-center timing)
  (define (act-count ov) (fxand act-overlap-mask ov))
  (define (pot-count ov) (fxand pot-overlap-mask ov))
  (let ([actthresh (fx* act-overlap-incr (tm-activation-threshold tm))]
        [redthresh (fx* act-overlap-incr (tm-reduced-basal-threshold tm))]
        [minthresh (tm-min-threshold tm)]
        [connected (tm-connected-permanence tm)]
        [iteration (tm-iteration tm)]
        [r-t-cells (list->fxvector reduced-threshold-cells)])
    (let next-source ([sources active-cells] [asegs (list)] [psegs (list)])
      ;; build active/matching segment lists [no benefit from building segx lists?]
      (if (pair? sources)
        (let* ( [source (car sources)]
                [segxv  (at-segxv at source)])
          (if segxv
            (let next-segment ([segxx (segxv-last segxv)] [asegs asegs] [psegs psegs])
              (if (fxpositive? segxx)
                (let* ( [segment (at-seg-ref at (segxv-ref segxv segxx))]
                        [synapse (synapses-search source (seg-synapses segment))])
                  (if synapse
                    (let ([new-overlap (fx1+ (fxmax iteration (seg-overlap segment)))]) ;; bump pot-count
                      (next-segment
                        (fx- segxx 2)
                        (if (fx<? (syn-perm synapse) connected)
                          (begin (seg-overlap-set! segment new-overlap)
                                 asegs)
                          (let ([new-overlap (fx+ new-overlap act-overlap-incr)])       ;; bump act-count
                            (seg-overlap-set! segment new-overlap)
                            (if (and (or (fx=? (act-count new-overlap) actthresh)
                                         (and (fx=? (act-count new-overlap) redthresh)
                                              (fxsearch r-t-cells (seg-cellx segment))))
                                     (not (memq segment asegs)))
                                (cons segment asegs)
                                asegs)))
                        (if (fx=? (pot-count new-overlap) minthresh)
                          (cons segment psegs)
                          psegs)))
                    (next-segment (fx- segxx 2) asegs psegs)))
                (next-source (cdr sources) asegs psegs)))
            (next-source (cdr sources) asegs psegs)))
        (values
          asegs
          (unique! eq?
            (sort-by! (lambda (sega segb)
                (fx<? (seg-cellx sega) (seg-cellx segb)))
              psegs)))))))
  ;(with-cost-center #t cost-center-1 csa-thunk))

;; === Supporting Functions ===
                                                                                            ;
(define (map-segments-to-cells segments) ;; {Seg} -> {CellX}
  ;; Produce sorted cellxs from segments
  (sort-unique! (map seg-cellx segments)))
                                                                                            ;
(define (map-sorted-segments-to-cells    ;; {Seg} -> {CellX}
          segments)
  (unique! fx=? (map seg-cellx segments)))
                                                                                            ;
(define (filter-segments-by-cell         ;; {Seg} {CellX} -> {Seg}
          segments cellxs)
  ;; Produce subset of segments that are on cellxs
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
(define (cols-from-cells tm cellxs)      ;; TM {CellX} -> {ColX}
  ;; Produce colxs of the cellxs (which are sorted), without duplicates
  (let ([cellx->colx (lambda (cellx) (cellx->colx tm cellx))])
    (unique! fx=? (map cellx->colx cellxs))))
                                                                                            ;
(define (cells-in-cols tm cellxs colxs)  ;; TM {CellX} {ColX} -> {CellX}
  ;; Produce cellxs for which (col of) cellx is in colxs
  (filter (lambda (cellx)
      (memv (cellx->colx tm cellx) colxs))
    cellxs))
                                                                                            ;
(define (exclude-segments tm segs colxs) ;; TM {Seg} {ColX} -> {Seg}
  ;; Produce segments for which cellx is not in colxs
  (filter (lambda (seg)
      (not (memv (cellx->colx tm (seg-cellx seg)) colxs)))
    segs))
                                                                                            ;
(define (filter-max-by-group tm          ;; TM {Seg} {CellX} (Seg -> Nat) -> {Seg}
          segs cells grouper)
  ;; Produce from segs, those with max potential overlaps in group, filtered by cells
  ;; (replaces filterSegmentsByCell/argmaxMulti/etc in .py)
  #| The JSP data structure diagram for segs is:  maybe-candidate-group*
                                                     /            \
                                              in-cells-groups°  not-in-cells-segs°
                                                    |              |
                                                  group*          seg*
                                                    |
                                                  seg*                          |#
  (let next-maybe-candidate-group ([segs segs] [results (list)])
    (cond
      [ (null? segs)  (reverse! results) ]
      [ else
        (let ([csegs (car segs)])
          (if (memv (seg-cellx csegs) cells)
            (let ([group (grouper csegs)])
              (let next-in-cells-seg ([segs (cdr segs)] [max-seg csegs]
                                      [max-overlap (seg-potential-overlaps tm csegs)])
                (cond
                  [ (null? segs)  (reverse! (cons max-seg results)) ]
                  [ (fx=? (grouper (car segs)) group)
                      (let ([overlap (seg-potential-overlaps tm (car segs))])
                        (if (fx>? overlap max-overlap)
                          (next-in-cells-seg (cdr segs) (car segs) overlap)
                          (next-in-cells-seg (cdr segs) max-seg max-overlap))) ]
                  [ else (next-maybe-candidate-group segs (cons max-seg results)) ])))
            (next-maybe-candidate-group (cdr segs) results))) ])))
                                                                                            ;
(define (all-cells-in-columns tm colxs)  ;; TM {ColX} -> {CellX}
  ;; Produce all cellxs in the colxs
  (let ([cpc (tm-cells-per-column tm)])
    (apply append! 
      (map (lambda (colx)
          (let ([first (fx* colx cpc)])
            (build-list cpc (lambda (i) (fx+ first i)))))
        colxs))))
                                                                                            ;
(define (get-axon-tree tm dend)          ;; TM Dendrite -> (Vectorof Source) (Vectorof {Seg})
  ;; Produce vectors of sources and segment lists of axon tree of given kind
  (let ([axon-tree (pba-axon-tree ((ba dend) tm))])
    (let-values ([(keys vals) (hashtable-entries (car axon-tree))])
      (values keys 
        (vector-map (lambda (val)
            (segxv-map (lambda (segx)
                (vector-ref (cdr axon-tree) segx))
              val))
          vals)))))
                                                                                            ;
(define (projection-of-source tm source  ;; TM Source PBA -> {ColX}
          pba)
  ;; Produce columns with cells which source synapses on; cache result for next use
  (let* ( [axon-tree (pba-axon-tree pba)]
          [segxv     (at-segxv axon-tree source)])
    (if segxv
      (let ([seg-table  (cdr axon-tree)]
            [latest-pos (cdr (hashtable-cell (pba-pos-cache pba) source (cons 2 (list))))]
            [lastx      (segxv-last segxv)])
        (do ( [segxx   (car latest-pos) (fx+ segxx 2)]
              [results (cdr latest-pos)
                  (let ([segment (vector-ref seg-table (segxv-ref segxv segxx))])
                    (cons (cellx->colx tm (seg-cellx segment)) results))])
            ((fx>? segxx lastx)
              (set-car! latest-pos segxx)
              (set-cdr! latest-pos results)
              results)))
      '())))
                                                                                            ;
(define (projection-of tm cellx pba)     ;; TM CellX PBA -> {ColX}
  ;; Produce columns with cells which cellx is connected to
  (projection-of-source tm (make-source (tm-cortical-column tm) (tm-layer tm) cellx) pba))

;; === Statistics ===
                                                                                            ;
(define (number-of-cells tm)             ;; TM -> Nat
  (fx* (tm-column-count tm) (tm-cells-per-column tm)))
                                                                                            ;
(define (number-of-connected-cells tm)   ;; TM -> Nat
  (number-of-cells tm))
                                                                                            ;
(define (number-of-segments tm dend)     ;; TM Dendrite -> Nat
  ;; produce number of segments (with more than 2 synapses)
  (let ([seg-index (cdr (pba-axon-tree ((ba dend) tm)))])
    (do ( [segx (vector-ref seg-index 0) (fx1- segx)]
          [nsegs 0 (fx+ nsegs
              (let ([seg (vector-ref seg-index segx)])
                (if seg  (if (fx>? (synapses-length (seg-synapses seg)) 2) 1 0)  0)))])
        ((fxzero? segx) nsegs))))
                                                                                            ;
(define (number-of-synapses tm dend)     ;; TM Dendrite -> Nat
  (let ([seg-index (cdr (pba-axon-tree ((ba dend) tm)))])
    (do ( [segx (vector-ref seg-index 0) (fx1- segx)]
          [nsyns 0 (fx+ nsyns
              (let ([seg (vector-ref seg-index segx)])
                (if seg  (synapses-length (seg-synapses seg))  0)))])
        ((fxzero? segx) nsyns))))
                                                                                            ;
(define (number-of-basal-segments tm)    ;; TM -> Nat
  (number-of-segments tm 'basal))
                                                                                            ;
(define (number-of-apical-segments tm)   ;; TM -> Nat
  (number-of-segments tm 'apical))
                                                                                            ;
(define (number-of-basal-synapses tm)    ;; TM -> Nat
  (number-of-synapses tm 'basal))
                                                                                            ;
(define (number-of-apical-synapses tm)   ;; TM -> Nat
  (number-of-synapses tm 'apical))
                                                                                            ;
(define (connection-lengths tm dend)     ;; TM Dendrite -> (Nat . Nat)
  ;; produce total length and number of connections of specified kind
  (let ([segxvs (hashtable-values (car (pba-axon-tree ((ba dend) tm))))])
    (let loop ([segxvx 0] [total 0] [count 0])
      (if (fx<? segxvx (vector-length segxvs))
        (let ([n (fxdiv (segxv-last (vector-ref segxvs segxvx)) 2)])
          (loop (fx1+ segxvx) (fx+ total n) (fx+ count (if (fxzero? n) 0 1))))
        (values total count)))))

)