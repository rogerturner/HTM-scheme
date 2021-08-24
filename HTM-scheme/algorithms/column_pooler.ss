;; === HTM-scheme Column Pooler algorithm  (C) 2019-21 Roger Turner. ===
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Based on code from Numenta Platform for Intelligent Computing (NuPIC) ;;
  ;; which is Copyright (C) 2016, Numenta, Inc.                            ;;
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
  
Translated from numenta htmresearch/.../column_pooler.py --
see comments there for descriptions of functions and parameters.
Follows logic of column_pooler.py except:
  cells in bursting minicolumns don't learn from feedforward input

*See htm_concept.ss for type and data structure descriptions and code conventions*
  
Selected comments from htmresearch/.../column_pooler.py:                                                                                       ;

This class constitutes a temporary implementation for a cross-column pooler.
The implementation goal of this class is to prove basic properties before
creating a cleaner implementation.

  maxSdrSize
    The maximum SDR size for learning.  If the column pooler has more
    than this many cells active, it will refuse to learn. This serves
    to stop the pooler from learning when it is uncertain of what object
    it is sensing.

  minSdrSize
    The minimum SDR size for learning. If the column pooler has fewer
    than this many active cells, it will create a new representation
    and learn that instead. This serves to create separate
    representations for different objects and sequences.
    If online learning is enabled, this parameter should be at least
    inertiaFactor*sdrSize.  Otherwise, two different objects may be
    incorrectly inferred to be the same, as SDRs may still be active
    enough to learn even after inertial decay.

  sampleSizeProximal
    Number of proximal synapses a cell should grow to each feedforward
    pattern, or -1 to connect to every active bit

  predictedInhibitionThreshold
    How much predicted input must be present for inhibitory behavior
    to be triggered. Only has effects if onlineLearning is true.

  sampleSizeDistal
    Number of distal synapses a cell should grow to each lateral
    pattern, or -1 to connect to every active bit

  inertiaFactor
    The proportion of previously active cells that remain
    active in the next timestep due to inertia (in the absence of
    inhibition). If onlineLearning is enabled, should be at most
    1 - learningTolerance, or representations may incorrectly become
    mixed.
    
  Learning mode: 

    If there was prior activity, we maintain it.

    If there are not enough previously active cells, then we are no longer on
    a familiar object. Either our representation decayed due to the passage
    of time (i.e. we moved somewhere else) or we were mistaken. Either way,
    create a new SDR and learn on it: randomly activate 'sdrSize' cells and create
    connections to incoming input.
    This case is the only way different object representations are created.
    These cells will represent the object and learn distal connections to each
    other and to lateral cortical columns.

    If we have a union of cells active, don't learn. This primarily affects
    online learning.

  Inference mode: 

    If there is some feedforward activity, perform
    spatial pooling on it to recognize previously known objects, then use
    lateral activity to activate a subset of the cells with feedforward
    support. 
    If there is no feedforward activity, use lateral activity to
    activate a subset of the previous active cells.

    First, activate the FF-supported cells that have the highest number of
    lateral active segments (as long as it's not 0)
    If we haven't filled the sdrSize quorum, add in inertial cells.

    We sort the previously-active cells by number of active lateral
    segments (this really helps). We then activate them in order of
    descending lateral activation.

    We use inertiaFactor to limit the number of previously-active cells
    which can become active, forcing decay even if we are below quota.

    Activate groups of previously active cells by order of their lateral
    support until we either meet quota or run out of cells.
    (on first touch activate all with feedforward support)

    If we haven't filled the sdrSize quorum, add cells that have feedforward
    support and no lateral support.

    Inhibit cells proportionally to the number of cells that have already
    been chosen. If ~0 have been chosen activate ~all of the feedforward
    supported cells. If ~sdrSize have been chosen, activate very few of
    the feedforward supported cells.

    Use the discrepancy:sdrSize ratio to determine the number of cells to
    activate.

  |#

  #!chezscheme

(library (HTM-scheme HTM-scheme algorithms column_pooler)
                                                                                            ;
(export
make-cp
compute
reset
get-active-cols
  cp-n-sdrs
  number-of-proximal-synapses
  number-of-proximal-segments
  number-of-distal-synapses
  number-of-distal-segments
  number-of-lateral-synapses
  number-of-lateral-segments
  (rename
    (cp-input-width   number-of-inputs)
    (cp-cell-count    number-of-cells)
    (cp-sdr-size      get-sdr-size)
    (cp-active-cells  get-active-cells))
  (rename
    (num-connected-proximal-synapses test:num-connected-proximal-synapses)
    (cp-min-sdr-size                 test:cp-min-sdr-size)
    (cp-max-sdr-size                 test:cp-max-sdr-size)
    (cp-proximal-permanences         test:cp-proximal-permanences)))
                                                                                            ;
(import
  (except (chezscheme) reset)
(HTM-scheme HTM-scheme algorithms htm_prelude)
(HTM-scheme HTM-scheme algorithms htm_concept)
(HTM-scheme HTM-scheme math       coordinates))
                                                                                            ;
  (implicit-exports #f)

;; === Layer record ===
                                                                                            ;
(define-record-type cp                   ;; CP
  (fields
    input-width
    lateral-input-widths
    cell-count
    sdr-size
    online-learning
    max-sdr-size
    min-sdr-size
    syn-perm-proximal-inc                ;; Proximal
    syn-perm-proximal-dec
    initial-proximal-permanence
    sample-size-proximal
    min-threshold-proximal
    connected-permanence-proximal
    predicted-inhibition-threshold
    syn-perm-distal-inc                  ;; Distal
    syn-perm-distal-dec
    initial-distal-permanence
    sample-size-distal
    activation-threshold-distal
    connected-permanence-distal
    inertia-factor                       ;; Fixnum3
    seed
    (mutable active-cells)
    proximal-permanences                 ;; CellVecOf Segment
    internal-distal-permanences          ;; CellVecOf Segment
    (mutable distal-permanences)         ;; {CellVecOf Segment}
    use-inertia
    cortical-column
    cells-per-column
    connect?
    (mutable n-sdrs))
  (sealed #t) (opaque #t) (nongenerative cp)
(protocol #;(make-cp kwargs)             ;; {KWarg} -> CP
  (lambda (new)
    (lambda (kwargs)
      (let* (
          [cp       (apply new (key-word-args kwargs cp-defaults))]
          [kwargs   (if (fxpositive? (cp-max-sdr-size cp))  kwargs
                        (cons `[max-sdr-size . ,(cp-sdr-size cp)] kwargs))]
          [kwargs   (if (fxpositive? (cp-min-sdr-size cp))  kwargs
                        (cons `[min-sdr-size . ,(cp-sdr-size cp)] kwargs))]
          [make-permanences
            (lambda ()
              (build-vector (cp-cell-count cp) (lambda (cellx)
                  (make-seg (cp-cortical-column cp) 0 cellx))))]
          [kwargs   (append `(
                        [proximal-permanences        . ,(make-permanences)]
                        [internal-distal-permanences . ,(make-permanences)]
                        [distal-permanences .
                          ,(build-list (length (cp-lateral-input-widths cp)) (lambda _
                                (make-permanences)))])
                      kwargs)])
        (random-seed (cp-seed cp))
        (apply new (key-word-args kwargs cp-defaults)))))))
                                                                                            ;
(define cp-defaults `(                   ;; (listof KWarg)
    [input-width                    . 0]
    [lateral-input-widths           . ()]
    [cell-count                     . 4096]
    [sdr-size                       . 40]
    [online-learning                . #f]
    [max-sdr-size                   . 0]
    [min-sdr-size                   . 0]
    [syn-perm-proximal-inc          . ,(perm 0.1)]
    [syn-perm-proximal-dec          . ,(perm 0.001)]
    [initial-proximal-permanence    . ,(perm 0.6)]
    [sample-size-proximal           . 20]
    [min-threshold-proximal         . 10]
    [connected-permanence-proximal  . ,(perm 0.50)]
    [predicted-inhibition-threshold . 20]
    [syn-perm-distal-inc            . ,(perm 0.1)]
    [syn-perm-distal-dec            . ,(perm 0.001)]
    [initial-distal-permanence      . ,(perm 0.6)]
    [sample-size-distal             . 20]
    [activation-threshold-distal    . 13]
    [connected-permanence-distal    . ,(perm 0.50)]
    [inertia-factor                 . ,(fx3<- 1.0)]
    [seed                           . 42]
    [active-cells                   . ()]
    [proximal-permanences           . #()]
    [internal-distal-permanences    . #()]
    [distal-permanences             . ()]
    [use-inertia                    . #t]
    [cortical-column                . 0]
    [cells-per-column               . 30]
    [connect?                       . #f]
    [n-sdrs                         . 0]))

;; === Column Pooler algorithm ===
                                                                                            ;
(define  compute                         ;; ...
  ;; run one time step of the column pooler algorithm
  (case-lambda
    [(cp feedforward-input learn)        ;; CP {CellX} Boolean ->
      (compute cp feedforward-input '() '() '() learn '() '()) ]
                                                                                            ;
    [(cp feedforward-input               ;; CP {CellX} {{CellX}} {CellX} Boolean ->
        lateral-inputs feedforward-growth-candidates learn)
      (compute cp feedforward-input '() lateral-inputs feedforward-growth-candidates learn '() '()) ]
                                                                                            ;
    [(cp feedforward-input distal-input  ;; CP {CellX} {{CellX}} {CellX} Boolean {CellX} ->
        lateral-inputs feedforward-growth-candidates learn predicted-input inhibited-cols)
      (let ([feedforward-growth-candidates
              (if (null? feedforward-growth-candidates) feedforward-input
                  feedforward-growth-candidates)])
        (if (cp-online-learning cp)
          (cond                          ;; online learning step
            [(fx>? (length predicted-input) (cp-predicted-inhibition-threshold cp))
              (let ([predicted-active-input
                      (intersect1d feedforward-input predicted-input)]
                    [predicted-growth-candidates
                      (intersect1d feedforward-growth-candidates predicted-input)])
                (compute-inference-mode cp predicted-active-input distal-input lateral-inputs)
                (compute-learning-mode cp predicted-active-input distal-input lateral-inputs
                  predicted-growth-candidates inhibited-cols)) ]  ; .py has feedforward-growth-candidates ?!
            [(not (fx<=? (cp-min-sdr-size cp)
                         (length (cp-active-cells cp))
                         (cp-max-sdr-size cp)))
              ;; if no single representation, try to infer one before attempting to learn
              (compute-inference-mode cp feedforward-input distal-input lateral-inputs)
              (compute-learning-mode cp feedforward-input distal-input lateral-inputs
                feedforward-growth-candidates inhibited-cols) ]
            [else 
              ;; no predicted input and single SDR, extend that representation
              (compute-learning-mode cp feedforward-input distal-input lateral-inputs
                feedforward-growth-candidates inhibited-cols) ])
          ;; not online-learning
          (if learn
                                         ;; learning step
            (compute-learning-mode cp feedforward-input distal-input lateral-inputs 
              feedforward-growth-candidates inhibited-cols)
                                         ;; inference step
            (compute-inference-mode cp feedforward-input distal-input lateral-inputs)))) ] ))
                                                                                            ;
(define (compute-learning-mode cp        ;; CP {CellX} {CellX} {{CellX}} {CellX} ->
          feedforward-input distal-input lateral-inputs feedforward-growth-candidates inhibited-cols)
  ;; in learning mode, maintain prior activity or create random sdr for new object
  (define (new-sdr)                      ;; -> SDR
    ;; produce random sdr with sdr-size bits
    (sort-unique-n! (u32-sample (iota (cp-cell-count cp)) (cp-sdr-size cp)) (cp-sdr-size cp)))
  ;; (compute-learning-mode)
  (let ([prev-active-cells (cp-active-cells cp)])
    (when (cond
        [ (fx<? (length prev-active-cells) (cp-min-sdr-size cp))
            ;; not enough previously active cells: create a new sdr and learn on it
            (cp-n-sdrs-set! cp (fx1+ (cp-n-sdrs cp)))
            (cp-active-cells-set! cp (new-sdr))
            #t ]
        [ (fx>? (length prev-active-cells) (cp-max-sdr-size cp))
            ;; union of cells active: don't learn
            #f ]
        [ else #t ] )
      ;; do the actual learning
      ;; (but [Behabadi etal 2012 Location-Dependent Excitatory Synaptic Interactions in
      ;;  Pyramidal Neuron Dendrites, Figure 8] reports that vertical L4 "driver" inputs
      ;;  terminate distally, horizontal [L2] connections proximally? )
      (learn! cp (cp-proximal-permanences cp)
        (setdiff1d prev-active-cells (cells-in-cols cp prev-active-cells inhibited-cols))
        feedforward-input feedforward-growth-candidates (cp-sample-size-proximal cp)
        (cp-initial-proximal-permanence cp) (cp-syn-perm-proximal-inc cp)
        (cp-syn-perm-proximal-dec cp))
      (do ( [i 0 (fx1+ i)]             ;; external distal learning
            [lateral-input lateral-inputs (cdr lateral-input)])
          ((fx=? i (length lateral-inputs)))
        (learn! cp (list-ref (cp-distal-permanences cp) i)
          (cp-active-cells cp) (car lateral-input) (car lateral-input)
          (cp-sample-size-distal cp) (cp-initial-distal-permanence cp)
          (cp-syn-perm-distal-inc cp) (cp-syn-perm-distal-dec cp)))
      (learn! cp (cp-internal-distal-permanences cp)
        (cp-active-cells cp) distal-input distal-input
        (cp-sample-size-distal cp) (cp-initial-distal-permanence cp)
        (cp-syn-perm-distal-inc cp) (cp-syn-perm-distal-dec cp)) )))
                                                                                            ;
(define (compute-inference-mode cp       ;; CP {CellX} {CellX} ->
          feedforward-input distal-input lateral-inputs)
  ;; in inference mode, if there is some feedforward activity, recognize previously
  ;; known objects and activate a subset of the cells with feedforward support;
  ;; otherwise use lateral activity to activate a subset of the previous active cells
  (define (extend-active active          ;; {CellX} {CellX} {Nat} Boolean -> {CellX}
            cells num-active-segs-for-cells allow-zero)
    ;; extend active with cells in order of activation until sdrSizequota is reached
    (let loop [ (ttop (apply fxmax num-active-segs-for-cells))
                (active active) ]
      (if (and (if allow-zero (not (fxnegative? ttop)) (fxpositive? ttop))
               (fx<? (length active) (cp-sdr-size cp)))
        (loop (fx1- ttop)
              (append
                (fold-left               ;; filter cells with activation = current ttop
                  (lambda (filtered cellx nasfc)
                    (if (fx=? nasfc ttop)
                        (cons cellx filtered)
                        filtered))
                  (list)
                  cells num-active-segs-for-cells)
                active))
        (list-sort fx<? active))))
  ;; (compute-inference-mode)
  (let* (
      [prev-active-cells (cp-active-cells cp)]
      [distal-input      (if (null? distal-input)  prev-active-cells
                             distal-input)]
      [overlaps          (compute-overlaps (cp-proximal-permanences cp)
                           feedforward-input (cp-connected-permanence-proximal cp))]
      [feedforward-supported-cells
                         (indices-where fx>=? overlaps (cp-min-threshold-proximal cp))]
      [num-active-segments-by-cell
                         (make-vector (cp-cell-count cp) 0)])
    (let (
        [overlaps (compute-overlaps (cp-internal-distal-permanences cp)
                      distal-input (cp-connected-permanence-distal cp))])
      (increment-where! num-active-segments-by-cell 
                        fx>=? overlaps (cp-activation-threshold-distal cp)))
    (when (pair? lateral-inputs)
      (for-each
        (lambda (lateral-input distal-permanence)
          (let ([overlaps (compute-overlaps distal-permanence
                            lateral-input (cp-connected-permanence-distal cp))])
            (increment-where! num-active-segments-by-cell 
                              fx>=? overlaps (cp-activation-threshold-distal cp))))
        lateral-inputs (cp-distal-permanences cp)))
    (let ([chosen-cells (list)])
      (unless (null? feedforward-supported-cells)
        (let ([num-active-segs-for-ff-sup-cells
                (list-of (vector-ref num-active-segments-by-cell x)
                         (x in feedforward-supported-cells))])
          (set! chosen-cells (extend-active chosen-cells  ;; exclude cells with 0 active segments
                               feedforward-supported-cells num-active-segs-for-ff-sup-cells #f))))
      (when (and (cp-use-inertia cp) (fx<? (length chosen-cells) (cp-sdr-size cp)))
        (let* ( [prev-cells   (setdiff1d prev-active-cells chosen-cells)]
                [inertial-cap (fxdiv (fx* (length prev-cells) (cp-inertia-factor cp)) fx3)])
          (when (fxpositive? inertial-cap)
            (let* ( [num-active-segs-for-prev-cells 
                      (list-of (vector-ref num-active-segments-by-cell x)
                               (x in prev-cells))]
                    [sort-indices (list-sort (lambda (x y)
                          (fx>? (list-ref num-active-segs-for-prev-cells x)
                                (list-ref num-active-segs-for-prev-cells y)))
                        (indexes num-active-segs-for-prev-cells))]
                    [prev-cells (take inertial-cap
                        (list-of (list-ref prev-cells x) (x in sort-indices)))]
                    [num-active-segs-for-prev-cells (take inertial-cap
                        (list-of (list-ref num-active-segs-for-prev-cells x) 
                                 (x in sort-indices)))])
              (set! chosen-cells (extend-active chosen-cells
                                   prev-cells num-active-segs-for-prev-cells #t))))))
      (let ([discrepancy (fx- (cp-sdr-size cp) (length chosen-cells))])
        (when (fxpositive? discrepancy)
          (let* ( [rem-ff-cells (setdiff1d feedforward-supported-cells chosen-cells)]
                  [n  (fxdiv (fx* (length rem-ff-cells) discrepancy) (cp-sdr-size cp))]
                  [n  (fxmin (fxmax n discrepancy) (length rem-ff-cells))])
            (if (fx>? (length rem-ff-cells) n)
              (let ([selected (u32-sample rem-ff-cells n)])
                (set! chosen-cells (append selected chosen-cells)))
              (set! chosen-cells (append rem-ff-cells chosen-cells))))))
      (cp-active-cells-set! cp (sort-unique! chosen-cells)))))
                                                                                            ;
(define (num-connected-proximal-synapses ;; CP {CellX} -> Nat
          cp cells)
  ;; produce total count of synapses above threshold of cells
  (let ([threshold (cp-connected-permanence-proximal cp)])
    (fold-left (lambda (total cellx)
        (fold-left (lambda (total synapse)
            (if (fx>=? (syn-perm synapse) threshold)
              (fx1+ total)
              total))
          total
          (seg-synapses->list (vector-ref (cp-proximal-permanences cp) cellx))))
      0
      cells)))
                                                                                            ;
(define  number-of-proximal-synapses     ;; CP {CellX} -> Nat
  ;; produce total count of synapses for cells
  (case-lambda 
    [(cp)
      (number-of-proximal-synapses cp (build-list (cp-cell-count cp) id))]
    [(cp cells)
      (fold-left (lambda (total cellx)
          (fx+ total
             (synapses-length (seg-synapses (vector-ref (cp-proximal-permanences cp) cellx)))))
        0
        cells)]))
                                                                                            ;
(define (number-of-proximal-segments cp) ;; CP -> Nat
  (vector-count (lambda (seg)
      (fxpositive? (synapses-length (seg-synapses seg))))
    (cp-proximal-permanences cp)))
                                                                                            ;
(define (number-of-distal-synapses cp)   ;; CP -> Nat
  (vector-fold-left (lambda (total seg)
      (fx+ total (synapses-length (seg-synapses seg))))
    0
    (cp-internal-distal-permanences cp)))
                                                                                            ;
(define (number-of-distal-segments cp)   ;; CP -> Nat
  (vector-count (lambda (seg)
      (fxpositive? (synapses-length (seg-synapses seg))))
    (cp-internal-distal-permanences cp)))
                                                                                            ;
(define (number-of-lateral-synapses cp)  ;; CP -> Nat
  (fold-left (lambda (total col)
      (fx+ total
        (vector-fold-left (lambda (total seg)
            (fx+ total (synapses-length (seg-synapses seg))))
          0
          col)))
    0
    (cp-distal-permanences cp)))
                                                                                            ;
(define (number-of-lateral-segments cp)  ;; CP -> Nat
  (fold-left (lambda (total col)
      (fx+ total
        (vector-count (lambda (seg)
            (fxpositive? (synapses-length (seg-synapses seg))))
          col)))
    0
    (cp-distal-permanences cp)))
                                                                                            ;
(define (reset cp)                       ;; CP ->
  ;; when learning this signifies we are to learn a unique new object
  (cp-active-cells-set! cp (list)))
                                                                                            ;
(define (cellx->colx cp cellx)           ;; CP CellX -> ColX
  (fxdiv cellx (cp-cells-per-column cp)))
                                                                                            ;
(define (cells-in-cols cp cellxs colxs)  ;; CP {CellX} {ColX} -> {CellX}
  ;; Produce cellxs for which (col of) cellx is in colxs
  (filter (lambda (cellx)
      (memv (cellx->colx cp cellx) colxs))
    cellxs))
                                                                                            ;
(define (get-active-cols cp)             ;; CP -> {ColX}
  ;; produce list of cols with active cell
  (unique! fx=? 
    (map (lambda (cellx)
        (cellx->colx cp cellx))
      (cp-active-cells cp))))
                                                                                            ;
(define (learn! cp permanences           ;; CP Permanences {CellX} {Source} {Source} Nat Perm Perm Perm ->
          active-cells active-input growth-candidate-input
          sample-size initial-permanence permanence-increment permanence-decrement)
  ;; for each active cell, reinforce active synapses, punish inactive synapses,
  ;; and grow new synapses to input bits that the cell isn't already connected to
  (let ([active-input-vec (list->fxvector active-input)])
    (adapt-synapses! permanences active-cells active-input-vec
                     permanence-increment permanence-decrement)
    (if (fxnegative? sample-size)
      (add-synapses!  cp permanences active-cells active-input initial-permanence)
      (grow-synapses! cp permanences active-cells active-input-vec growth-candidate-input
                      sample-size initial-permanence))))

;; === Synapses and Permanences ===
;; ('Permanences' in column_pooler.py correspond to Segments in HTM-scheme)
                                                                                            ;
(define (adapt-synapses! permanences     ;; Permanences {CellX} {Source} Perm Perm ->
          active-cells active-input-vec permanence-increment permanence-decrement)
  ;; update synapses: strengthen those connected to input, weaken others, remove on zero
  (for-each (lambda (cellx)
      (let ([segment (vector-ref permanences cellx)])
        (seg-synapses-set! segment
          (synapses-update! (lambda (synapse)
              (let* ( [source  (syn-source synapse)]
                      [perm    (if (fxsearch active-input-vec source)
                                 (clip-max (fx+ (syn-perm synapse) permanence-increment))
                                 (fx- (syn-perm synapse) permanence-decrement))])
                (and (fxpositive? perm) (make-syn source perm))))
            (seg-synapses segment)))))
    active-cells))
                                                                                            ;
(define (add-synapses! cp permanences    ;; CP Permanences {CellX} {Source} Perm ->
          active-cells active-input initial-perm)
  ;; create new synapses on all active-cells for all connectable input bits
  (for-each (lambda (active-cell)  
      (let* (
          [post-mcx         (cellx->colx cp active-cell)]
          [connect?         (cp-connect? cp)]
          [sources
            (if connect?
              (filter (lambda (source)
                  (connect? source post-mcx))
                active-input)
              active-input)]
          [segment (vector-ref permanences active-cell)])
        (seg-synapses-set! segment
          (synapses-add!
            sources
            (length sources)
            initial-perm
            (seg-synapses segment)))))
    active-cells))
                                                                                            ;
(define (grow-synapses! cp permanences   ;; CP Permanences {CellX} (FXVectorof Source) {Source} {Nat} Perm ->
          active-cells active-input-vec growth-candidates sample-size initial-perm)
  ;; create new synapses for some growth-candidates ('setRandomZerosOnOuter' in .py)
  (for-each (lambda (active-cell)
      (let* (
          [segment   (vector-ref permanences active-cell)]
          [synapses  (seg-synapses segment)]
          [max-new
            (do ( [sx    (fx1- (synapses-length synapses)) (fx1- sx)]
                  [count 0 (if (fxsearch active-input-vec (syn-source (synapses-ref synapses sx)))
                              (fx1+ count)
                              count)] )
                ((fxnegative? sx) (fx- sample-size count)))])
        (when (fxpositive? max-new)
          (let* (
              [connect?  (cp-connect? cp)]
              [post-mcx  (cellx->colx cp active-cell)]
              [sources
                (if connect?
                  (filter (lambda (source)
                      (connect? source post-mcx))
                    growth-candidates)
                  growth-candidates)]
              [n-new   (length sources)])
            (when (fxpositive? n-new)
              (seg-synapses-set! segment
                (synapses-add!
                  (if (fx<=? n-new max-new)  sources
                      (sort-unique-n! (u32-sample sources max-new) max-new))
                  (fxmin n-new max-new)
                  initial-perm
                  synapses)))))))
    active-cells))
                                                                                            ;
(define (compute-overlaps                ;; Permanences {InputX} Perm -> (CellVecOf Nat)
          permanences inputxs threshold)
  ;; produce counts of overlaps between synapses above threshold and inputs
  ;; ('rightVecSumAtNZGteThresholdSparse' in htmresearch/column_pooler.py)
  (vector-map (lambda (segment)
      (let ([synapses (seg-synapses segment)])
        (fold-left (lambda (overlaps inputx)
            (let ([synapse (synapses-search inputx synapses)])
              (if (and synapse (fx>=? (syn-perm synapse) threshold))
                (fx1+ overlaps)
                overlaps)))
          0
          inputxs)))
    permanences))
    
;; === Vector and list utilities ===
                                                                                            ;
(define (indices-where pred? vec value)  ;; (X Y -> Bool) (Vectorof X) Y -> {Nat}
  ;; produce sorted list of indices of vec elements for which (pred? element value)
  (do ( [vx (fx1- (vector-length vec)) (fx1- vx)]
        [result (list) (if (pred? (vector-ref vec vx) value)
                           (cons vx result)
                           result)] )
      ((fxnegative? vx) result)))
                                                                                            ;
(define (increment-where! cs pred? xs y) ;; (Vectorof Nat) (X Y -> Boolean) (Vectorof X) Y ->
  ;; mutates cs, incrementing elements corresponding to xs for which (pred? x y)
  (do ( [vx 0 (fx1+ vx)] )
      ((fx=? vx (vector-length xs)))
      (when (pred? (vector-ref xs vx) y)
        (vector-set! cs vx (fx1+ (vector-ref cs vx))))))
                                                                                            ;
(define-syntax fold-of                   ;; Op Base Expr Clause ...
  ;; List comprehensions from the Programming Praxis standard prelude
  ;; (see https://programmingpraxis.com/contents/standard-prelude/ for details)
  ;; Only (list-of (var in list)) is used
  (syntax-rules (range in is)
    ((_ "z" f b e) (set! b (f b e)))
    ((_ "z" f b e (v range fst pst stp) c ...)
      (let* ((x fst) (p pst) (s stp)
             (le? (if (positive? s) <= >=)))
        (do ((v x (+ v s))) ((le? p v) b)
          (fold-of "z" f b e c ...))))
    ((_ "z" f b e (v range fst pst) c ...)
      (let* ((x fst) (p pst) (s (if (< x p) 1 -1)))
        (fold-of "z" f b e (v range x p s) c ...)))
    ((_ "z" f b e (v range pst) c ...)
      (fold-of "z" f b e (v range 0 pst) c ...))
    ((_ "z" f b e (x in xs) c ...)
      (do ((t xs (cdr t))) ((null? t) b)
        (let ((x (car t)))
          (fold-of "z" f b e c ...))))
    ((_ "z" f b e (x is y) c ...)
      (let ((x y)) (fold-of "z" f b e c ...)))
    ((_ "z" f b e p? c ...)
      (if p? (fold-of "z" f b e c ...)))
    ((_ f i e c ...)
      (let ((b i)) (fold-of "z" f b e c ...)))))
                                                                                            ;
(define-syntax list-of (syntax-rules ()  ;; Expr Clause ...
  ((_ arg ...) (reverse (fold-of
    (lambda (d a) (cons a d)) '() arg ...)))))
                                                                                            ;
(define-syntax sum-of (syntax-rules ()   ;; Expr Clause ...
  ((_ arg ...) (fold-of + 0 arg ...))))

)