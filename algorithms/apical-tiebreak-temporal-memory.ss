#| HTM-scheme Apical Tiebreak Temporal Memory (C) 2019-2022 Roger Turner.
   License: AGPL3 https://www.gnu.org/licenses/agpl-3.0.txt (see Notices below)

This HTM-scheme implementation of the Apical Tiebreak Temporal Memory algorithm
aims to replicate key behaviour of Numenta implementations[1] using a different
underlying data structure:
• Numenta:    SparseMatrixConnections (C++, Python/numpy: see tutorial) [2]
• HTM-scheme: packed 32-bit synapses, Scheme Hashtables for connections [3]

 [1] github.com/numenta/nupic.research/tree/master/packages/columns,
     ~/src/nupic/research/frameworks/columns/apical_tiebreak_temporal_memory.py
 [2] ~nupic.research.core/blob/master/examples/bindings/sparse_matrix_how_to.py
     ~nupic.research.core/blob/master/src/nupic/math/SparseMatrix.hpp
 [3] github.com/rogerturner/HTM-scheme/blob/master/frameworks/htm-concept.ss
     (see introductory comments for an overview of HTM-scheme data structures)

For ease of comparison with the Numenta algorithm ("attm.py") Scheme code uses
matching function and variable names and is structured similarly. (For example 
"_calculateBasalLearning" in attm.py calculates "learningActiveBasalSegments",
"calculate-basal-learning" below calculates "learning-active-basal-segments".)
See comments in attm.py for further description of functions/parameters below.

Indentation facilitates using an editor "Fold All" view for a file overview.

Description from attm.py:
  A generalized Temporal Memory with apical dendrites that add a "tiebreak".
  Basal connections are used to implement traditional Temporal Memory.
  The apical connections are used for further disambiguation. If multiple cells
  in a minicolumn have active basal segments, each of those cells is predicted,
  unless one of them also has an active apical segment, in which case only the
  cells with active basal and apical segments are predicted.

Follows logic of attm.py except:
  activate-cells has an optional bursting-columns input (see layer4.ss for usage)
  a connect? hook in grow-synapses can be used to filter growth-candidates
  grow-synapses can prefer growth candidates already connected near to segment
  destroy-threshold determines when synapses are removed from segments
  winnerCells in attm.py renamed learning-cells

"Remember that all models are wrong; the practical question is how wrong
do they have to be to not be useful." [George Box]
  
  |#

  #!chezscheme

(library (algorithms apical-tiebreak-temporal-memory)
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
  tm-cells-per-column
  at-basalht
  at-apicalht
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
  cells-in-cols
  cols-from-cells
  )
                                                                                            ;
(import
  (except (chezscheme) reset)
          (frameworks htm-prelude)
          (frameworks htm-concept))
                                                                                            ;
  (implicit-exports #f)
  
;; === Types (see htm-concept.ss for common HTM types) ===
                                                                                            ;
;; Dendrite     = Symbol: 'proximal | 'basal | 'apical
;; TM           = Record: tm parameters
;; PBA          = Record: proximal/basal/apical parameters

;; === Parameter defaults, data structures, accessors ===
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
    [destroy-threshold                     . ,(perm 0.0)]  ;; default value => attm.py behaviour
    [use-bursting-columns-input            . #f]           ;; ''
    [basal-connect?                        . #f]           ;; ''
    [apical-connect?                       . #f]           ;; ''
    [prefer-nearby-minicol                 . #f] ))        ;; ''
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
    seed                                 ;; Integer
    (mutable active-cells)               ;; (Listof CellX)
    (mutable learning-cells)             ;; (Listof CellX)
    (mutable predicted-cells)            ;; (Listof CellX)
    (mutable predicted-active-cells)     ;; (Listof CellX)
    proximal-ds                          ;; PBA [proximal dendrite parameters]
    basal-ds                             ;; PBA [basal dendrite parameters]
    apical-ds                            ;; PBA [apical dendrite parameters]
    use-apical-tiebreak                  ;; T?
    use-apical-modulation-basal-threshold;; T?
    ;; following parameters are not in apical_tiebreak_temporal_memory.py
    (mutable iteration)                  ;; Nat [incrementing by #x10000]
    cortical-column                      ;; Nat
    layer                                ;; Nat
    axon-radius2                         ;; Nat
    destroy-threshold                    ;; Perm
    use-bursting-columns-input           ;; T?
    basal-connect?                       ;; (CellX ColX -> T?) | #f
    apical-connect?                      ;; (CellX ColX -> T?) | #f
    prefer-nearby-minicol )              ;; T?
  (opaque #t) (nongenerative tm)
(protocol #;(make-tm kwargs)             ;; {KWarg} -> TM
  ;; use kwargs to override defaults, construct initial values, create record
  (lambda (new)
    (lambda (kwargs)
      (let* (
          [tm       (apply new (key-word-args kwargs attm-defaults))]
          [n-cells  (fx* (tm-column-count tm) (tm-cells-per-column tm))]
          [kwargs   (cons `[basal-ds . ,(make-pba
                              (tm-basal-input-size tm)
                              (tm-basal-predicted-segment-decrement tm)
                              (make-at (fxdiv n-cells 4))
                              (tm-basal-connect? tm)
                              (make-bytevector n-cells 0)
                              (list)
                              (list))] kwargs)]
          [kwargs   (cons `[apical-ds . ,(make-pba
                              (tm-apical-input-size tm)
                              (tm-apical-predicted-segment-decrement tm)
                              (make-at (fxdiv n-cells 4))
                              (tm-apical-connect? tm)
                              #f
                              (list)
                              (list))] kwargs)])
        (apply new (key-word-args kwargs attm-defaults)))))))
                                                                                            ;
(define-record-type pba (fields          ;; PBA [proximal/basal/apical parameters]
    ;; initialised in make-tm
    input-size                           ;; Nat
    predicted-segment-decrement          ;; Perm
    axon-tree                            ;; AxonTree
    connect?                             ;; (Source ColX -> T?) | #f
    segment-count                        ;; (Bytevector CellX->U8) | #f
    (mutable active-segments)            ;; {Segment}
    (mutable matching-segments))         ;; {Segment}
  (sealed #t) (opaque #t) (nongenerative pba))
                                                                                            ;
(define (at-basalht tm)                  ;; TM -> Hashtable
  ;; accessor for basal axon-tree hashtable
  (car (pba-axon-tree (tm-basal-ds tm))))
                                                                                            ;
(define (at-apicalht tm)                 ;; TM -> Hashtable
  ;; accessor for apical axon-tree hashtable
  (car (pba-axon-tree (tm-apical-ds tm))))
                                                                                            ;
(define (seg-count+ pba cellx)           ;; PBA Cellx ->
  ;; increment segment count for cellx
  (let ([sc (pba-segment-count pba)])
    (when sc
      (bytevector-u8-set! sc cellx (fxmin 255 (fx1+ (bytevector-u8-ref sc cellx)))))))
                                                                                            ;
(define (seg-count- pba cellx)           ;; PBA Cellx ->
  ;; decrement segment count for cellx
  (let ([sc (pba-segment-count pba)])
    (when sc
      (bytevector-u8-set! sc cellx (fxmax 0 (fx1- (bytevector-u8-ref sc cellx)))))))
                                                                                            ;
(define (seg-count@ pba cellx)           ;; PBA Cellx -> Nat
  ;; produce segment count for cellx (0 if no count maintained)
  (let ([sc (pba-segment-count pba)])
    (if sc  (bytevector-u8-ref sc cellx)  0 )))
                                                                                            ;
(define (ba dend)                        ;; Dendrite -> (TM -> PBA)
  ;; produce basal-ds or apical-ds 
  (case dend
    [(basal) tm-basal-ds]
    [else    tm-apical-ds]))
                                                                                            ;
(define-syntax with-tm
  ;; "(define (fn tm) / (with-tm / (form)))" defines cellx->colx in form using fn's tm
  (lambda (x)
    (syntax-case x ()
    [(k e1 e2 ...)
      (with-syntax ([tm          (datum->syntax #'k 'tm)]
                    [cellx->colx (datum->syntax #'k 'cellx->colx)])
        #'(let ([cpc (tm-cells-per-column tm)])
            (define (cellx->colx cellx) (fxdiv cellx cpc))
            e1 e2 ... )) ])))

;; === Apical Tiebreak Temporal Memory Algorithm ===
                                                                                            ;
(define (reset tm)                       ;; TM ->
  ;; clear all cell and segment activity (overlap counts are stored in segments)
  (tm-active-cells-set!           tm (list))
  (tm-learning-cells-set!         tm (list))  ;; htmresearch "winner" renamed learning
  (tm-predicted-cells-set!        tm (list))
  (tm-predicted-active-cells-set! tm (list))
  (pba-active-segments-set!   (tm-basal-ds  tm) (list))
  (pba-active-segments-set!   (tm-apical-ds tm) (list))
  (pba-matching-segments-set! (tm-basal-ds  tm) (list))
  (pba-matching-segments-set! (tm-apical-ds tm) (list)))
                                                                                            ;
(define (depolarize-cells tm             ;; TM {ColX} {CellX} {CellX} T? ->
          active-columns                 ;; minicolumns with proximal input
          basal-input                    ;; active cells pre-synaptic to basal segments
          apical-input                   ;; active cells pre-synaptic to apical segments
          learning)                      ;; #f (inference mode) to use reduced-basal-threshold
  ;; save predicted cells, active and matching segments (overlaps saved in segments)
  (tm-iteration-set! tm (fx+ iteration-incr (tm-iteration tm)))
                                         ;; (see calculate-segment-activity for iteration-incr use)
  (when (fxpositive? (tm-cells-per-column tm))
    (let*-values (
        [(active-cols-mask) (target-cols->mask active-columns)]
        [(active-apical-segments matching-apical-segments)
            (calculate-segment-activity tm apical-input '() (tm-apical-ds tm) active-cols-mask)]
        [(reduced-basal-threshold-cells)
           (if (or learning (not (tm-use-apical-modulation-basal-threshold tm)))  '()
               (map-segments-to-cells active-apical-segments))]
        [(active-basal-segments matching-basal-segments)
            (calculate-segment-activity tm basal-input reduced-basal-threshold-cells
                                          (tm-basal-ds tm) active-cols-mask)])
      (tm-predicted-cells-set! tm
        (calculate-predicted-cells tm active-basal-segments active-apical-segments))
      (pba-active-segments-set!   (tm-basal-ds  tm) active-basal-segments)
      (pba-active-segments-set!   (tm-apical-ds tm) active-apical-segments)
      (pba-matching-segments-set! (tm-basal-ds  tm) matching-basal-segments)
      (pba-matching-segments-set! (tm-apical-ds tm) matching-apical-segments))))
                                                                                            ;
(define (activate-cells tm               ;; TM {ColX} {CellX} {CellX} {CellX} {CellX} T? . {ColX} ->
          active-columns                 ;; minicolumns with proximal input
          basal-reinforce-candidates     ;; pre-synaptic cells whose basal synapses should be incremented
          apical-reinforce-candidates    ;; pre-synaptic cells whose apical synapses should be incremented
          basal-growth-candidates        ;; potential pre-synaptic cells for basal segments
          apical-growth-candidates       ;; potential pre-synaptic cells for apical segments
          learning . bursting-columns)   ;; update segments/synapses? optional bursting columns input
  ;; activate cells in the specified columns, using predictions from depolarize-cells, then learn.
  ;;   active-cells:           correct-predicted-cells + all cells in bursting columns
  ;;   learning-cells:         correct-predicted-cells + cells of learning-matching-basal-segments + new-basal-segment-cells
  ;;   predicted-active-cells: correct-predicted-cells
  ;; update synapses on active and matching segments, and create new segments
(define (learn pba                       ;; PBA {Seg} {Seg} {Seg} {CellX} {Source} {Source} ->
          learning-active-segments learning-matching-segments
          segments-to-punish new-segment-cells
          reinforce-candidates growth-candidates)
  ;; adapt learning segments, punish incorrect predictions, grow new segments
  (let ([learning-segments    (append! learning-active-segments learning-matching-segments)]
        [reinforce-candidates (list->u32-vector reinforce-candidates)]
        [increment            (tm-permanence-increment tm)]
        [decrement            (tm-permanence-decrement tm)]
        [destroy              (tm-destroy-threshold tm)])
    (define (update synapse)             ;; Synapse -> Synapse | #f
      ;; produce synapse with adjusted permanence, or #f if below threshold
      (let ([source (syn-source synapse)]
            [perm   (syn-perm synapse)])
        (if (u32-search reinforce-candidates source)
          (make-syn source (clip-max (fx+ perm increment)))
          (let ([perm (fx- perm decrement)])
            (and (fx>? perm destroy) (make-syn source perm))))))
    ;; (learn)                           ;; learn on existing segments
    (adapt-and-grow-segments tm learning-segments growth-candidates pba update)
    (let ([decrement (fx- (pba-predicted-segment-decrement pba))])
      (unless (fxzero? decrement)        ;; punish incorrect predictions
        (punish-segments tm segments-to-punish reinforce-candidates decrement pba)))
    (when (pair? growth-candidates)
      (let* (                            ;; grow new segments
          [at (pba-axon-tree pba)]
          [new-segments (mapr (lambda (cellx)
                            (seg-count+ pba cellx)
                            (at-make-seg at cellx))
                          new-segment-cells) ])
        (adapt-and-grow-segments tm new-segments growth-candidates pba #f)))))
  ;; (activate-cells)                    ;; calculate active cells
  (when (fxpositive? (tm-cells-per-column tm))
    (let* (
        [active-cols-mask (target-cols->mask active-columns)]
        [correct-predicted-cells
          (cells-in-active-cols tm (tm-predicted-cells tm) active-cols-mask)]
        [bursting-columns
          (if (and (tm-use-bursting-columns-input tm) (pair? bursting-columns))
            (car bursting-columns)       ;; bursting-columns, if supplied, is ({ColX} . '())
            (setdiff1d active-columns (cols-from-cells tm (tm-predicted-cells tm))))]
        [new-active-cells
          (append correct-predicted-cells (all-cells-in-columns tm bursting-columns))])
      (if learning
        (let*-values (                   ;; calculate learning
            [ ( learning-active-basal-segments
                learning-matching-basal-segments
                basal-segments-to-punish
                new-basal-segment-cells
                learning-cells)
                  (calculate-basal-learning tm
                    active-cols-mask bursting-columns correct-predicted-cells) ]
            [ ( learning-active-apical-segments
                learning-matching-apical-segments
                apical-segments-to-punish
                new-apical-segment-cells)
                  (calculate-apical-learning tm learning-cells active-cols-mask) ] )
          (learn (tm-basal-ds tm) learning-active-basal-segments learning-matching-basal-segments
            basal-segments-to-punish new-basal-segment-cells
            basal-reinforce-candidates basal-growth-candidates)
          (learn (tm-apical-ds tm) learning-active-apical-segments learning-matching-apical-segments
            apical-segments-to-punish new-apical-segment-cells
            apical-reinforce-candidates apical-growth-candidates)
                                         ;; save the results
          (tm-learning-cells-set! tm learning-cells))
        (tm-learning-cells-set! tm '())) ;; (not learning)
      (tm-active-cells-set!           tm (sort-unique! new-active-cells))
      (tm-predicted-active-cells-set! tm correct-predicted-cells))))
                                                                                            ;
(define (calculate-basal-learning tm     ;; TM {ColX} {ColX} {CellX} -> {Seg} {Seg} {Seg} {CellX} {CellX} 
          active-cols-mask bursting-columns correct-predicted-cells)
  ;; basic Temporal Memory learning: see comments in apical_tiebreak_temporal_memory.py
  (let* (
      [learning-active-basal-segments
        (filter-segments-by-cell! (pba-active-segments (tm-basal-ds tm)) correct-predicted-cells)]
      [matching-basal-segments
        (pba-matching-segments (tm-basal-ds tm))]
      [cells-for-matching-basal
        (map-segments-to-cells matching-basal-segments)]
      [matching-cells-in-bursting-columns
        (cells-in-cols tm cells-for-matching-basal bursting-columns)]
      [bursting-columns-with-no-match
        (setdiff1d bursting-columns (cols-from-cells tm cells-for-matching-basal))]
      [learning-matching-basal-segments
        (choose-best-segment-per-column tm matching-basal-segments matching-cells-in-bursting-columns)]
      [new-basal-segment-cells
        (get-cells-with-fewest-segments tm bursting-columns-with-no-match)]
      [learning-cells (sort-unique!
          (append correct-predicted-cells new-basal-segment-cells
                  (map seg-cellx learning-matching-basal-segments)))]
      [basal-segments-to-punish
        (segments-not-active! tm matching-basal-segments active-cols-mask)])
    (values
      learning-active-basal-segments
      learning-matching-basal-segments
      basal-segments-to-punish
      new-basal-segment-cells
      learning-cells)))
                                                                                            ;
(define (calculate-apical-learning tm    ;; TM {CellX} {ColX} -> {Seg} {Seg} {Seg} {CellX}
          learning-cells active-cols-mask)
  ;; calculate apical learning for each learning cell - see comments in apical_tiebreak_temporal_memory.py
  (let* (
      [learning-active-apical-segments
        (filter-segments-by-cell! (pba-active-segments (tm-apical-ds tm)) learning-cells)]
      [learning-cells-without-active-apical
        (setdiff1d learning-cells (map-segments-to-cells learning-active-apical-segments))]
      [matching-apical-segments
        (pba-matching-segments (tm-apical-ds tm))]
      [cells-for-matching-apical
        (map-segments-to-cells matching-apical-segments)]
      [learning-cells-with-matching-apical
        (intersect1d learning-cells-without-active-apical cells-for-matching-apical)]
      [learning-matching-apical-segments
        (choose-best-segment-per-cell tm learning-cells-with-matching-apical matching-apical-segments)]
      [new-apical-segment-cells
        (setdiff1d learning-cells-without-active-apical learning-cells-with-matching-apical)]
      [apical-segments-to-punish
        (segments-not-active! tm matching-apical-segments active-cols-mask)])
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
  ;; calculate predicted cells for given active segments
  (with-tm
    (let ([cells-for-basal-segments (map-segments-to-cells active-basal-segments)])
      (if (not (tm-use-apical-tiebreak tm))
        cells-for-basal-segments
        (let* (
            [cells-for-apical-segments (map-segments-to-cells active-apical-segments)]
            [fully-depolarized-cells   (intersect1d cells-for-basal-segments
                                                    cells-for-apical-segments)]
            [partly-depolarized-cells  (setdiff1d cells-for-basal-segments
                                                  fully-depolarized-cells)]
            [inhibited-mask            (in1d (map cellx->colx partly-depolarized-cells)
                                             (map cellx->colx fully-depolarized-cells))])
          (union1d fully-depolarized-cells
                    (exclude-by-mask partly-depolarized-cells inhibited-mask)))))))
                                                                                            ;
#;(learn, learn-on-new-segments)         ;; see learn in activate-cells
                                                                                            ;
(define (choose-best-segment-per-cell tm ;; TM {CellX} {Seg} -> {Seg}
          cells all-matching-segments)
  ;; choose matching segment with max active potential synapses; previous code (cf attm.py):
  ;; (let* ([candidate-segments  (filter-segments-by-cell all-matching-segments cells)]
  ;;        [one-per-cell-filter (argmax-multi (map (lambda (s) (potential-overlaps tm s)) candidate-segments)
  ;;                                           (map seg-cellx candidate-segments))])
  ;;        [learning-segments   (list-refs! candidate-segments one-per-cell-filter)])
  (filter-max-by-group tm all-matching-segments cells seg-cellx))
                                                                                            ;
(define (choose-best-segment-per-column  ;; TM {Seg} {CellX} -> {Seg}
          tm all-matching-segments matching-cells)
  ;; like choose-best-segment-per-cell but for cells in col; previous code (cf attm.py):
  ;; (let* ([candidate-segments     (filter-segments-by-cell all-matching-segments matching-cells)]
  ;;        [cell-scores            (map (lambda (s) (potential-overlaps tm s)) candidate-segments)]
  ;;        [columns-for-candidates (map (lambda (s) (cellx->colx tm (seg-cellx s))) candidate-segments)]
  ;;        [one-per-column-filter  (argmax-multi cell-scores columns-for-candidates)])
  ;;      (list-refs! candidate-segments one-per-column-filter))
  (with-tm
    (filter-max-by-group tm all-matching-segments matching-cells
      (lambda (seg) (cellx->colx (seg-cellx seg))))))
                                                                                            ;
(define (get-cells-with-fewest-segments  ;; TM {ColX} -> {CellX}
          tm colxs)
  ;; produce cellx for each colx with fewest basal segments; break ties randomly
  (map (lambda (colx)
      (let* ( [start (fx* colx (tm-cells-per-column tm))]
              [limit (fx+ start (tm-cells-per-column tm))])
        (let loop ([cellx start] [cellxs (list)] [fewest (greatest-fixnum)])
          (if (fx<? cellx limit)
            (let ([n-segs (seg-count@ (tm-basal-ds tm) cellx)]
                  [nextx  (fx1+ cellx)])
              (cond
                [ (fx<? n-segs fewest) (loop nextx (list cellx)        n-segs) ]
                [ (fx=? n-segs fewest) (loop nextx (cons cellx cellxs) n-segs) ]
                [ else                 (loop nextx cellxs              fewest) ] ))
            (if (null? (cdr cellxs))  (car cellxs)
                (list-ref cellxs (random (length cellxs))))))))
    colxs))
                                                                                            ;
(define (random-cell-in-columns tm colxs);; TM {ColX} -> {CellX}
  ;; produce random cell in each specified column (alternative to get-cells-with-fewest-segments)
  (let ([cpc (tm-cells-per-column tm)])
    (map (lambda (colx)
        (let ([first (fx* colx cpc)])
          (fx+ first (random cpc))))
      colxs)))

;; === Accessors (getActiveCells etc) === (provided by define-record-type)

;; === Connections ===
                                                                                            ;
(define (punish-segments tm segments     ;; TM {Segment} (U32VectorOf Source) Perm PBA ->
          reinforce-candidates permanence-delta pba)
  ;; update permanences of synapses connected to reinforce-candidates
  ;; (replaces connections.adjustActiveSynapses)
(define (remove-reference source seg)    ;; Source Segment -> #f
  ;; remove source's reference to seg from axon tree
  (let* ( [axon-tree (pba-axon-tree pba)]
          [target    (at-target axon-tree source)])
    (when target
      (segxv-remove! (lambda (segx)
          (eq? seg (at-seg-ref axon-tree segx)))
        target)))
  #f )
  ;; (punish-segments)
  (let ([destroy-threshold (tm-destroy-threshold tm)])
    (let next-seg ([segs segments])
      (when (pair? segs)
        (let ([segment (car segs)])
          (synapses-update! (lambda (synapse)
              (let ([source (syn-source synapse)])
                (if (u32-search reinforce-candidates source)
                  (let ([perm (fx+ (syn-perm synapse) permanence-delta)])
                    ;; (used to punish incorrect predictions, so negative delta only)
                    (if (fx<=? perm destroy-threshold)
                      (remove-reference source segment)  ;; => #f
                      (make-syn source perm)))
                  synapse)))
            segment)
          (when (fx<? (synapses-length (seg-synapses segment)) (tm-min-threshold tm))
            (synapses-for-each (lambda (s)
                (remove-reference (syn-source s) segment))
              segment)
            (at-free-seg (pba-axon-tree pba) (seg-segx segment))
            (seg-count- pba (seg-cellx segment))))
        (next-seg (cdr segs))))))
                                                                                            ;
(define (adapt-and-grow-segments tm      ;; TM {Segment} {Source} PBA (Synapse -> Synapse) ->
          segments growth-candidates pba update)
  ;; add synapses for selection of growth-candidates to segments, update existing synapses
  ;; (replaces connections.adjustSynapses and connections.growSynapsesToSample)
  (define (grow-synapses segment         ;; Segment ColX {Source} {Source} ->
             colx preferred others)
    ;; add synapses to segment (which is in colx) that could connect
    (define (connectable colx sources)   ;; ColX {Source} -> {Source}
      ;; exclude sources by connect rules, or already connected
      (let ([connect? (pba-connect? pba)])
        (if connect?
          (filter (lambda (gc)           ;; filter => connect? is called on gcs reversed by pairs
              (if (fxnegative? colx)
                (and (connect? gc colx) (not (synapses-search gc segment)))
                (connect? gc colx)))
            sources)
          (remp (lambda (gc) (synapses-search gc segment)) sources))))
    ;; (grow-synapses):
    (let* (                              ;; prefer-nearby-minicol #f => preferred is '()
        [preferred    (connectable -1 preferred)]
        [n-preferred  (length preferred)] 
        [max-new      (if (fxnegative? (tm-sample-size tm))
                        (greatest-fixnum)  ;; (potential-overlaps is zero for new segments)
                        (fxmax 0 (fx- (tm-sample-size tm)
                                      (potential-overlaps tm segment))))]
        [max-new      (if (fxnegative? (tm-max-synapses-per-segment tm))  max-new
                          (fxmin max-new
                                 (fx- (tm-max-synapses-per-segment tm)
                                      (synapses-length (seg-synapses segment)))))]
        [preferred    (if (fx<=? n-preferred max-new)  preferred
                              (list-sample preferred max-new n-preferred))]
        [max-new      (fx- max-new n-preferred)]
        [others       (let* ( [others (connectable colx others)]
                              [n-others (length others)])
                        (if (fx<=? n-others max-new)  others
                            (list-sample others max-new n-others)))]
        [sources      (union1d! preferred others)])
      (at-update (pba-axon-tree pba) colx sources (tm-initial-permanence tm) segment update)))
  ;; (adapt-and-grow-segments):
  (with-tm
    (for-each (lambda (seg)
        (let ([colx (cellx->colx (seg-cellx seg))])
          (let-values ([(preferred others)
                        (if (tm-prefer-nearby-minicol tm)
                          (!partition (lambda (gc)
                              (let ([target (at-target (pba-axon-tree pba) gc)])
                                (and target (target-adjacent? target colx))))
                            growth-candidates)
                          (values '() growth-candidates))])
            (grow-synapses seg colx preferred others))))
      segments)))
                                                                                            ;
  ;; pack overlap counts and iteration "time-stamp" into one segment field
  (define pot-overlap-mask #x000FF)
  (define act-overlap-mask #x0FF00)
  (define act-overlap-incr #x00100)
  (define iteration-incr   #x10000)
                                                                                            ;
(define (potential-overlaps tm seg)      ;; TM Segment -> Nat
  ;; produce potential overlaps for seg (ignoring a value from a previous iteration)
  (let ([overlap (seg-overlap seg)])
    (if (fx<? overlap (tm-iteration tm))  0
        (fxand pot-overlap-mask overlap))))
                                                                                            ;
(define (calculate-segment-activity tm   ;; TM {Source} {CellX} PBA -> {Seg} {Seg}
          active-input reduced-basal-threshold-cells pba active-cols-mask)
  ;; for each source in active-input use basal/apical axon tree to find segments with
  ;; synapses for that input, build segment lists and update overlap counts per thresholds
  ;; (combines attm.py calculateApical/BasalSegmentActivity and connections.computeActivity)
  (define (act-count ov) (fxand act-overlap-mask ov))
  (define (pot-count ov) (fxand pot-overlap-mask ov))
  (let* ( [at        (pba-axon-tree pba)]
          [actthresh (fx* act-overlap-incr (tm-activation-threshold tm))]
          [redthresh (fx* act-overlap-incr (tm-reduced-basal-threshold tm))]
          [minthresh (tm-min-threshold tm)]
          [connected (tm-connected-permanence tm)]
          [iteration (tm-iteration tm)]
          [seg-table (cdr at)])
    (let next-source ([sources active-input] [asegs (list)] [psegs (list)])
      ;; build active/matching segment lists [no benefit from building segx lists?]
      (if (pair? sources)
        (let* ( [source (car sources)]
                [target (at-target at source)])
          ;; only active segments are used to determine learning, which is only on active cols
          ;; (cells with proximal input), so only sources targeting active cols need => active segments
          (if (and target (target-and? target active-cols-mask))
            (let ([segxx-last (segxv-last target)])
              (let next-segment ([segxx segxv-base] [asegs asegs] [psegs psegs])
                (if (fx<=? segxx segxx-last)
                  (let* ( [segment (vector-ref seg-table (segxv-ref target segxx))]
                          [synapse (synapses-search source segment)])
                    (if synapse
                      (let ([new-overlap  ;; overlap with bumped pot-count, for current iteration
                              (fx1+ (fxmax iteration (seg-overlap segment)))])
                        (next-segment
                          (fx+ segxx segx-bytes)
                          (if (fx>=? (syn-perm synapse) connected)
                            (let ([new-overlap (fx+ new-overlap act-overlap-incr)])
                              (seg-overlap-set! segment new-overlap)  ;; save overlap with new counts
                              (let ([new-act (act-count new-overlap)])
                                (cond
                                  [(fx>? new-act actthresh)  asegs ]
                                  [(and (fx=? new-act actthresh)
                                        (not (seg-memq segment asegs)))
                                      (cons segment asegs) ]
                                  [(and (fx=? new-act redthresh)
                                        (sdr-memv (seg-cellx segment) reduced-basal-threshold-cells)
                                        (not (seg-memq segment asegs)))
                                      (cons segment asegs) ]
                                  [else asegs ])))
                            (begin (seg-overlap-set! segment new-overlap)  ;; save overlap (non-connected)
                              asegs))
                          (if (fx=? (pot-count new-overlap) minthresh)
                            (cons segment psegs)
                            psegs)))
                      (next-segment (fx+ segxx segx-bytes) asegs psegs)))
                  (next-source (cdr sources) asegs psegs))))
            ;; matching segments only
            (if target
              (let ([segxx-last (segxv-last target)])
                (let next-segment ([segxx segxv-base] [psegs psegs])
                  (if (fx<=? segxx segxx-last)
                    (let* ( [segment (vector-ref seg-table (segxv-ref target segxx))]
                            [new-overlap (fx1+ (fxmax iteration (seg-overlap segment)))])
                      (seg-overlap-set! segment new-overlap)             ;; save overlap with new pot-count
                      (next-segment
                        (fx+ segxx segx-bytes)
                        (if (fx=? (pot-count new-overlap) minthresh)
                          (cons segment psegs)
                          psegs)))
                    (next-source (cdr sources) asegs psegs))))
              (next-source (cdr sources) asegs psegs))))
          (values
            asegs
            (sort-unique-by! (lambda (sega segb)
                (let ([cellxa (seg-cellx sega)] [cellxb (seg-cellx segb)])
                  (cond                    ;; canonical order of segs for uniquing
                    [(fx<? cellxa cellxb) ]
                    [(fx=? cellxa cellxb)
                       (fx<? (seg-segx sega) (seg-segx segb)) ]
                    [else #f ])))
              psegs))))))

;; === Supporting Functions ===
                                                                                            ;
(define (sdr-memv x xs)                  ;; Fixnum {Fixnum} -> {Fixnum} | #f
  ;; produce first tail of xs (which is sorted) with car equal to x, or #f
  (let next-x ([xs xs])
    (cond
      [(null? xs)  #f ]
      [(fx<? (car xs) x) (next-x (cdr xs)) ]
      [(fx=? (car xs) x) xs ]
      [else #f ])))
                                                                                            ;
(define (seg-memq seg segs)              ;; Segment {Segment} -> {Segment} | #f
  ;; produce first tail of segs with car eq? seg, or #f
  (let next-seg ([segs segs])
    (cond
      [(null? segs) #f ]
      [(eq? seg (car segs)) segs ]
      [else (next-seg (cdr segs)) ])))
                                                                                            ;
(define (mapr proc xs)                   ;; (X -> Y) {X} -> {Y}
  ;; produce proc applied to each x (result reversed)
  (let next-x ([xs xs] [ys (list)])
    (cond
      [(null? xs) ys ]
      [else (next-x (cdr xs) (cons (proc (car xs)) ys)) ])))
                                                                                            ;
(define (map-segments-to-cells segments) ;; {Seg} -> {CellX}
  ;; produce sorted cellxs from segments
  (sort-unique! (mapr seg-cellx segments)))
                                                                                            ;
(define (filter-segments-by-cell!        ;; {Seg}! {CellX} -> {Seg}
          segs cellxs)
  ;; reduce segments to those that are on cellxs
  (filter! (lambda (seg)
      (sdr-memv (seg-cellx seg) cellxs))
    segs))
                                                                                            ;
(define (map-segments-to-synapse-counts  ;; {Seg} -> {Nat}
          segs)
  (mapr (lambda (seg)
      (synapses-length (seg-synapses seg)))
    segs))
                                                                                            ;
(define (cellx->colx tm cellx)           ;; TM CellX -> ColX
  (fxdiv cellx (tm-cells-per-column tm)))
                                                                                            ;
(define (cols-from-cells tm cellxs)      ;; TM {CellX} -> {ColX}
  ;; produce sorted colxs of the cellxs (which are sorted), without duplicates
  (with-tm
    (sort-unique! (map cellx->colx cellxs))))
                                                                                            ;
(define (cells-in-cols tm cellxs colxs)  ;; TM {CellX} {ColX} -> {CellX}
  ;; produce cellxs for which (col of) cellx is in colxs
  (with-tm
    (filter (lambda (cellx)
        (sdr-memv (cellx->colx cellx) colxs))
      cellxs)))
                                                                                            ;
(define (cells-in-active-cols tm cellxs  ;; TM {CellX} ColsMask -> {CellX}
          active-cols-mask)
  ;; produce cellxs for which (col of) cellx is in active-cols
  (with-tm
    (filter (lambda (cellx)
        (bytevector-bit-set? active-cols-mask (cellx->colx cellx)))
      cellxs)))
                                                                                            ;
(define (segments-not-active! tm segs    ;; TM {Seg}! ColsMask -> {Seg}
          active-cols-mask)
  ;; reduce segs to those for which cellx is not in active cols
  (with-tm
    (filter! (lambda (seg)
        (not (bytevector-bit-set? active-cols-mask (cellx->colx (seg-cellx seg)))))
      segs)))
                                                                                            ;
(define (filter-max-by-group tm          ;; TM {Seg} {CellX} (Seg -> Nat) -> {Seg}
          segs cells grouper)
  ;; produce from segs, those with max potential overlaps in group, filtered by cells
  ;; (replaces filterSegmentsByCell/argmaxMulti/etc in attm.py)
  #| The Jackson structure diagram for segs is:   maybe-candidate-group*
                                                     /            \
                                            in-cells-groupsº  not-in-cells-segsº
                                                    |              |
                                                  group*          seg*
                                                    |
                                                   seg*                          |#
  (let next-maybe-candidate-group ([segs segs] [results (list)])
    (cond
      [ (null? segs)  (!reverse! results) ]
      [ else
        (let ([csegs (car segs)])
          (if (sdr-memv (seg-cellx csegs) cells)
            (let ([group (grouper csegs)])
              (let next-in-cells-seg ([segs (cdr segs)] [max-seg csegs]
                                      [max-overlap (potential-overlaps tm csegs)])
                (cond
                  [ (null? segs)  (!reverse! (cons max-seg results)) ]
                  [ (fx=? (grouper (car segs)) group)
                      (let ([overlap (potential-overlaps tm (car segs))])
                        (if (fx>? overlap max-overlap)
                          (next-in-cells-seg (cdr segs) (car segs) overlap)
                          (next-in-cells-seg (cdr segs) max-seg max-overlap))) ]
                  [ else (next-maybe-candidate-group segs (cons max-seg results)) ])))
            (next-maybe-candidate-group (cdr segs) results))) ])))
                                                                                            ;
(define (all-cells-in-columns tm colxs)  ;; TM {ColX} -> {CellX}
  ;; produce all cellxs in the colxs
  (let ([cpc (tm-cells-per-column tm)])
    (apply append! 
      (map (lambda (colx)
          (let ([first (fx* colx cpc)])
            (build-list cpc (lambda (i) (fx+ first i)))))
        colxs))))
                                                                                            ;
(define (get-axon-tree tm dend)          ;; TM Dendrite -> (Vectorof Source) (Vectorof {Seg})
  ;; produce vectors of sources and segment lists of axon tree of given kind
  (let ([axon-tree (pba-axon-tree ((ba dend) tm))])
    (let-values ([(keys vals) (hashtable-entries (car axon-tree))])
      (values keys 
        (vector-map (lambda (val)
            (segxv-map (lambda (segx)
                (vector-ref (cdr axon-tree) segx))
              val))
          vals)))))

;; === Statistics ===
                                                                                            ;
(define (number-of-cells tm)             ;; TM -> Nat
  (fx* (tm-column-count tm) (tm-cells-per-column tm)))
                                                                                            ;
(define (number-of-connected-cells tm)   ;; TM -> Nat
  (number-of-cells tm))
                                                                                            ;
(define (number-of-segments tm dend)     ;; TM Dendrite -> Nat
  (let ([seg-table (cdr (pba-axon-tree ((ba dend) tm)))])
    (do ( [segx (fx1- (vector-length seg-table)) (fx1- segx)]
          [nsegs 0 (fx+ nsegs
              (let ([seg (vector-ref seg-table segx)])
                (if (fixnum? seg)  0  1)))])
        ((fxzero? segx) nsegs))))
                                                                                            ;
(define (number-of-synapses tm dend)     ;; TM Dendrite -> Nat
  (let ([seg-table (cdr (pba-axon-tree ((ba dend) tm)))])
    (do ( [segx (fx1- (vector-length seg-table)) (fx1- segx)]
          [nsyns 0 (fx+ nsyns
              (let ([seg (vector-ref seg-table segx)])
                (if (fixnum? seg)  0  (synapses-length (seg-synapses seg)))))])
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
(define (connection-lengths tm dend)     ;; TM Dendrite -> Nat Nat
  ;; produce total length and number of connections of specified kind
  (let ([segxvs (hashtable-values (car (pba-axon-tree ((ba dend) tm))))])
    (let loop ([segxvx 0] [total 0] [count 0])
      (if (fx<? segxvx (vector-length segxvs))
        (let ([n (fxdiv (segxv-last (vector-ref segxvs segxvx)) segx-bytes)])
          (loop (fx1+ segxvx) (fx+ total n) (fx+ count (if (fxzero? n) 0 1))))
        (values total count)))))

)

#| Notices:

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    
  Contact: https://discourse.numenta.org/u/rogert   |#
