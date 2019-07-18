#!chezscheme

;; === HTM-scheme Untangling Sequences project Copyright 2019 Roger Turner. ===
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

  ;; Partial replication and extension of experiments reported in Numenta paper
  ;; Ahmad & Hawkins 2017 "Untangling Sequences: Behavior vs. External Causes"
  ;; doi: 10.1101/190678 
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.

(import
          (except (chezscheme) add1 make-list random reset)
          (HTM-scheme HTM-scheme algorithms htm_prelude)
          (HTM-scheme HTM-scheme algorithms htm_concept)
  (prefix (HTM-scheme HTM-scheme algorithms L2objL4locL4seq_patch)         l2l4tm:)
  (prefix (HTM-scheme HTM-scheme algorithms apical_tiebreak_temporal_memory) attm:)
  (prefix (HTM-scheme HTM-scheme algorithms column_pooler)                     cp:))

;; === Types (see htm_concept.ss) ===
;; SDR            = Listof Nat
;; CColX          = Range num-cortical-columns
;; FeatureSDRX    = Range num-features
;; LocationSDRX   = Range num-locations
;; Sensation      = (FeatureSDRX, LocationSDRX|#f)
;; Object         = Vectorof Sensation
;; Sequence       = Vectorof (FeatureSDRX, #f)
;; Experience     = Object|Sequence

(define-record-type sensation (fields
  feature                                ;; SDRX
  location))                             ;; SDRX

(define (experiment exp-args run-args)   ;; KWargs KWargs ->
  ;; train l2l4tm-patch on object/sequence experiences, test, summarize response
  (let* [
    (args (append run-args exp-args))
    (get  (lambda (key default)
            (let ((specified (assoc key args)))
              (if specified (cdr specified) default))))
                                                                                            ;
;; model and learning parameter defaults (overridable by experiment/run parameters)
    (num-cortical-columns         (get 'num-cortical-columns          1))
    (column-count                 (get 'column-count                512))
    (num-cells-per-column         (get 'num-cells-per-column         16))
    (num-input-bits               (get 'num-input-bits (min 20 (max 4 
                                          (int<- (* column-count 0.04))))))
    (initial-permanence           (get 'initial-permanence   (perm 0.41)))
    (connected-permanence         (get 'connected-permanence (perm 0.6 )))
    (permanence-increment         (get 'permanence-increment (perm 0.1 )))
    (permanence-decrement         (get 'permanence-decrement (perm 0.02)))
    (activation-threshold         (get 'activation-threshold (min (- num-input-bits 1)
                                          (int<- (* num-input-bits 0.9)))))
    (min-threshold                (get 'min-threshold (min (- num-input-bits 1)
                                          (int<- (* num-input-bits 0.9)))))
    (reduced-basal-threshold      (get 'reduced-basal-threshold
                                           (int<- (* activation-threshold 0.6))))
    (predicted-segment-decrement  (get 'predicted-segment-decrement (perm 0.001)))
    (enable-feedback              (get 'enable-feedback      #t))
    (online-learning              (get 'online-learning      #f))
                                                                                            ;
;; setup network: override algorithm default parameters (scaled to values above)
    (L2-overrides  (get 'L2-overrides `(
      [cell-count                         . ,(min 4096 
                                              (int<- (/ (* column-count num-cells-per-column) 2)))]
      [sample-size-proximal               . ,(int<- (* 0.75 num-input-bits))]
      [min-threshold-proximal             . ,(int<- (* 0.5  num-input-bits))]
      [predicted-inhibition-threshold     . ,(int<- (* 0.8  num-input-bits))]
      [sample-size-distal                 . ,num-input-bits]
      [activation-threshold-distal        . ,(int<- (* 0.6  num-input-bits))]
      [sdr-size                           . ,(* 2 num-input-bits)]
      [syn-perm-proximal-dec              . ,(perm 0.001)]
      [online-learning                    . ,online-learning])))
    (TM-overrides  (get 'TM-overrides `(
      [initial-permanence                 . ,initial-permanence]
      [connected-permanence               . ,connected-permanence]
      [permanence-increment               . ,permanence-increment]
      [permanence-decrement               . ,permanence-decrement  #;(perm 0.03) ]
      [activation-threshold               . ,activation-threshold]
      [min-threshold                      . ,min-threshold]
      [reduced-basal-threshold            . ,reduced-basal-threshold]
      [apical-predicted-segment-decrement . ,predicted-segment-decrement  #;(perm 0.0) ]
      [basal-predicted-segment-decrement  . ,predicted-segment-decrement])))
    (L4-overrides  (get 'L4-overrides `(
      [initial-permanence                 . ,initial-permanence  #;(perm 0.21) ]
      [connected-permanence               . ,connected-permanence]
      [permanence-increment               . ,permanence-increment]
      [permanence-decrement               . ,permanence-decrement]
      [activation-threshold               . ,activation-threshold]
      [min-threshold                      . ,min-threshold]
      [reduced-basal-threshold            . ,reduced-basal-threshold]
      [apical-predicted-segment-decrement . ,predicted-segment-decrement  #;(perm 0.0) ]
      [basal-predicted-segment-decrement  . ,predicted-segment-decrement])))
    (patch 
      (l2l4tm:make-patch num-cortical-columns column-count num-input-bits
        enable-feedback L2-overrides L4-overrides TM-overrides))
                                                                                            ;
;; default experiment parameters (overridable)
    (figure                            (get 'figure               'z))
    (num-objects                       (get 'num-objects          50))
    (num-points                        (get 'num-points           10))
    (num-sequences                     (get 'num-sequences        50))
    (seq-length                        (get 'seq-length           10))
    (num-features                      (get 'num-features        100))
    (num-locations                     (get 'num-locations       100))
    (num-repetitions                   (get 'num-repetitions       5))
    (settling-time                     (get 'settling-time         1))
    (random-seq-location               (get 'random-seq-location  #f))
    (no-seq-location                   (get 'no-seq-location      #f))
    (intersperse-noise                 (get 'intersperse-noise     0))
    (superimpose-sequence              (get 'superimpose-sequence #f))
    (interleave-training               (get 'interleave-training  #f))
                                                                                            ;
;; create the sequences and objects
    (random-sdr 
      (lambda (size bits)         ;; Nat -> (X -> SDR)
        ;; produce function to generate SDR of size with sparsity bits/size
        (let ((range (build-vector size id)))
          (lambda _ (list-sort fx<? (vector->list 
              (vector-sample range bits)))))))
    (external-input-size
      (if (<= 512 column-count 1024) 1024
          (* 2 column-count)))
    (feature-pool                      ;; Vectorof SDR indexed by FeatureSDRX
      (build-vector num-features  (random-sdr column-count num-input-bits)))
    (location-pool                     ;; Vectorof SDR indexed by LocationSDRX
      (build-vector num-locations (random-sdr external-input-size num-input-bits)))
    (random-locations                  ;; Vectorof SDR indexed by LocationSDRX
      (build-vector num-locations (if no-seq-location
                                      (lambda _ (list))
                                      (random-sdr external-input-size num-input-bits))))
    (create-random-experiences         ;; Nat Nat Nat -> Experiences
      (lambda (num-exps num-sensations num-locations)
      ;; produce Object experiences, or Sequence experiences when num-locations is zero
      (let ((location-array (if (zero? num-locations) #f
                (build-vector (max num-locations num-sensations)
                  (if (>= num-locations num-sensations)
                    id
                    (lambda _ (random num-locations)))))))
        (build-vector num-exps (lambda _
          ;; locations are distinct but features can be repeated
          (let ((locations (if (positive? num-locations)
                               (vector-sample location-array num-sensations)
                               (make-vector num-sensations #f))))
            (build-vector num-sensations (lambda (sx)
                (make-sensation (random num-features) (vector-ref locations sx))))))))))
    (sequences                         ;; Vectorof Experience
      (create-random-experiences num-sequences seq-length 0))
    (objects                           ;; Vectorof Experience
      (create-random-experiences num-objects num-points num-locations))
                                                                                          ;
    ;; Stats is (ExpVectorOf (SensVectorOf (KeyListOf (CCVectorOf Number|{Nat}))))
    (stats          (build-vector (max num-objects num-sequences)
                      (lambda _ 
                        (build-vector (if (eq? figure 'f6)  (* 8 num-points)
                                          (max num-points seq-length))
                                      (lambda _ (list))))))
    (sequence-order (apply append (make-list num-repetitions (build-list num-sequences id))))
    (sequence-order (vector-sample (list->vector sequence-order) (* num-repetitions num-sequences)))
    (display-timing (get 'display-timing  #t))
    (train-keys     (get 'train-keys     '()))
    (test-keys      (get 'test-keys      '()))
    (start-time     (cpu-time))
    ]

#;> (define (have-sensations-of          ;; Experience Bool ( -> (SDR Nat -> SDR)) [(Nat -> )] ->
              experience learn seq . report)
      ;; feed each sensation of experience settling-time times, report after each sensation
      (let* ( 
          (object?    (sensation-location (vector-ref experience 0)))
          (ns         (vector-length experience))
          (sensations (if object?        ;; shuffle object sensations
                          (vector-sample experience ns)
                          experience))
          ;; for each sensation of object, ccs get different selection of feature+location
          ;; for sequences, all ccs get same feature
          (selects    (build-vector ns (lambda (sx)
                          (if object?
                            (build-vector num-cortical-columns (lambda (ccx)
                                (modulo (+ sx ccx) ns)))
                            (make-vector num-cortical-columns sx))))))
        (do ((sx 0 (add1 sx))) ((= sx ns))
          (let* ( 
              (select (vector-ref selects sx))
              (features
                (build-vector num-cortical-columns (lambda (ccx)
                  (if object?
                    (let ((object-sdr (vector-ref feature-pool
                            (sensation-feature (vector-ref sensations (vector-ref select ccx))))))
                      (if (and superimpose-sequence seq)
                        (unique! fx=? (union1d object-sdr 
                                 (vector-ref feature-pool 
                                    (sensation-feature (vector-ref (vector-ref sequences seq) sx)))))
                        object-sdr))
                    (vector-ref feature-pool 
                                (sensation-feature (vector-ref sensations sx)))))))
              (locations
                (if object?
                  (vector-refs location-pool
                    (build-vector num-cortical-columns (lambda (ccx)
                      (sensation-location (vector-ref sensations (vector-ref select ccx))))))
                  (build-vector num-cortical-columns (lambda _
                      (vector-ref (if (or no-seq-location random-seq-location)
                                      random-locations 
                                      location-pool) 
                                  (random num-locations)))))))
            (when object?                ;; settling-time applies to L4 location population only
              (do ((_ 0 (add1 _))) ((>= _ (- settling-time 1)))
                (l2l4tm:compute patch features locations (or learn online-learning) (l2l4tm:L4 loc))))
            (l2l4tm:compute patch features locations (or learn online-learning) (l2l4tm:L4 loc seq))
            (unless (null? report)
              ((car report) sx features))))))
                                                                                            ;
#;> (define (get-stats-for keys)         ;; {Symbol} [{{Nat}}]-> (Alistof (CCVectorOf {Nat}))
      ;; produce alist of stats for keys
      (define (lengths ->cells ->layers) ;; (Layer -> {CellX}) (Patch -> CCVectorOf Layer) -> (vectorof Nat)
        ;; produce lengths of cells of layer for each cc
        (vector-map length
          (vector-map ->cells (->layers patch))))
      (map (lambda (key)
          (cons key
            (case key
              [(L4lnp) (lengths attm:get-predicted-cells        l2l4tm:patch-L4s)]
              [(L4lpa) (lengths attm:get-predicted-active-cells l2l4tm:patch-L4s)]
              [(L4ac)  (vector-map attm:get-active-cells           (l2l4tm:patch-L4s patch))]
              [(L4pa)  (vector-map attm:get-predicted-active-cells (l2l4tm:patch-L4s patch))]
              [(TMlnp) (lengths attm:get-predicted-cells        l2l4tm:patch-TMs)]
              [(TMlpa) (lengths attm:get-predicted-active-cells l2l4tm:patch-TMs)]
              [(TMac)  (vector-map attm:get-active-cells           (l2l4tm:patch-TMs patch))]
              [(TMpa)  (vector-map attm:get-predicted-active-cells (l2l4tm:patch-TMs patch))]
              )))
        keys))
                                                                                            ;
#;> (define (train-all es)               ;; Experiences ->
      ;; train on each e num-repetitions times
      (do-with-progress (vector-length es)
        (lambda (ex)
          (let ((e (vector-ref es ex)))
            (do ((r 0 (add1 r))) ((>= r (- num-repetitions 1)))
              (let ((seq (if superimpose-sequence
                            (vector-ref sequence-order (+ (* ex num-repetitions) r))
                            #f)))
                (have-sensations-of e #t seq)
                                           ;; reset TM between presentations of sequence
                (unless (sensation-location (vector-ref e 0))
                  (vector-for-each attm:reset (l2l4tm:patch-TMs patch)))))
            (let ((seq (if superimpose-sequence
                          (vector-ref sequence-order (+ (* ex num-repetitions) (- num-repetitions 1)))
                          #f)))
              (have-sensations-of e #t seq ;; last training repetition: save stats
                (lambda (sx features)
                  (vector-set! (vector-ref stats ex) sx 
                    (append (get-stats-for train-keys) (list (cons 'feat features)) '())))))
            (when (positive? intersperse-noise)
              (l2l4tm:reset patch))
            (let ((noise (vector-ref (create-random-experiences 1 seq-length 0) 0)))
              (do ((_ 0 (add1 _))) ((= _ intersperse-noise))
                (have-sensations-of noise #t #f)))
            (l2l4tm:reset patch)))))     ;; reset after each object
                                                                                            ;
#;> (define (train-all-interleaved)      ;; ->
      ;; train objects & sequences, interleaving in (sqrt num-repetitions) batches
      ;; reset only after inner batch of each experience
      (assert (= num-objects num-sequences))
      (let*-values (
        [(n-inner remainder) (exact-integer-sqrt num-repetitions)]
        [(n-outer)           (if (>= remainder n-inner)
                                 (add1 n-inner)
                                 n-inner)])
      (do-with-progress n-outer
        (lambda (o)
          (vector-for-each
            (lambda (ex)
              (do ((i 0 (add1 i))) ((fx= i n-inner))
                (let ((seq (if superimpose-sequence
                              (vector-ref sequence-order (+ (* ex num-repetitions)
                                                            (+ (* o n-inner) i)))
                              #f)))
                  (unless superimpose-sequence
                    (have-sensations-of (vector-ref sequences ex) #t #f))
                  (have-sensations-of (vector-ref objects   ex) #t seq)))
              (l2l4tm:reset patch))
            (indexes objects))))))
                                                                                            ;
#;> (define (infer-all es)               ;; Experiences ->
      ;; test each experience and save stats for each sensation
      (do-with-progress (vector-length es)
        (lambda (ex)
          (have-sensations-of (vector-ref es ex) #f #f
            (lambda (sx features)
              (vector-set! (vector-ref stats ex) sx
                (append (get-stats-for test-keys)
                        (vector-ref (vector-ref stats ex) sx)
                        '()))))
          (l2l4tm:reset patch))))
                                                                                            ;
#;> (define (infer-switching)            ;; ->
      ;; test mix of objects & sequences
      (let ((objectxs (vector-sample (build-vector num-objects id) 8)))
        (for-each
          (lambda (item-type ix)
            (have-sensations-of
              (vector-ref
                (if (eq? item-type 'seq) sequences  objects) (vector-ref objectxs ix))
              #f #f
              (lambda (sx features)
                (vector-set! (vector-ref stats 0) (+ (* 10 ix) sx)
                  (get-stats-for test-keys)))))
          '(seq obj seq obj seq seq obj seq) (build-list 8 id))))
                                                                                            ;
;; run experiment: 1 train network
    (if interleave-training
        (train-all-interleaved)
      (begin
        (unless superimpose-sequence
          (train-all sequences))
        (train-all objects)))
;;                 2 reset and run inference
    (l2l4tm:reset patch)
    (case figure
      [(f6)
        (infer-switching) ]
      [else
        (infer-all objects)
        (infer-all sequences) ] )
;;                 3 save statistics for plots
    (let ((summary-file "HTM-scheme/projects/untangling_sequences/experiment.data")
          (mean  (lambda (x) (if (eq? figure 'f6)  x
                         (int<- (/ x (+ num-objects num-sequences)))))))
      (with-output-to-file summary-file (lambda ()
        (write [cons* 
          [list "figure" (symbol->string figure)]
          [list "using"  run-args]
          (map (lambda (key)             ;; each stats key (Symbol -> )
            [cons (symbol->string key)
              ;; stats is (ExpVectorOf (SensVectorOf (KeyListOf (CCVectorOf Number))))
              [list (vector->list (vector-map mean   ;; average of exps (for 4a/5a)
                (vector-fold-left 
                  (lambda (accs sens-for-exp)  ;; (listof Number) (SensListOf (KeyListOf (vectorof Number)))
                    ;; ..accumulated objects|sequences
                    (if (null? (vector-ref sens-for-exp 0))  accs
                      (vector-map (lambda (ac keys-stats)    ;; ..each sensation, average ccs
                          ;; Number (Alistof (vectorof Number)) -> Number
                        (let ((stats-for-key (assq key keys-stats)))
                          (define (add-overlap-with train-key)
                            (let ((trained (assq train-key keys-stats)))
                              (if trained
                                (+ ac (vector-average 
                                        (vector-map (lambda (t a)
                                            (length (unique! fx=? (intersect1d t a))))
                                          (cdr trained)
                                          (cdr stats-for-key))))
                                ac)))
                          (if stats-for-key
                            (case key
                              [(L4pa)  (add-overlap-with 'L4ac)]
                              [(TMpa)  (add-overlap-with 'TMac)]
                              [(TMlnp TMlpa L4lnp L4lpa)
                                (+ ac (vector-average (cdr stats-for-key)))]
                              [else ac])
                            ac)))
                        accs 
                        sens-for-exp)))        ;; -> (vectorof Number)
                  (build-vector (vector-length (vector-ref stats 0))
                    (lambda _ 0))
                  stats)) ) ] ] )
            test-keys) ] ))
            'truncate ))
      (when display-timing
        (let ((tenths (fxdiv (+ (- (cpu-time) start-time) 50) 100))
              (sum-ccs (lambda (f tms)
                (vector-fold-left (lambda (sum tm)
                    (+ sum (f tm)))
                  0
                  (tms patch)))))
          (for-each display `(
              ,(symbol->string figure) " " ,@run-args #\newline
              ,(fxdiv tenths 10) "." ,(fxmod tenths 10) "s thread time\n"
              ,(sum-ccs attm:tm-n-segments-created l2l4tm:patch-L4s) " L4 segments\n"
              ,(sum-ccs attm:tm-n-synapses-created l2l4tm:patch-L4s) " L4 synapses\n"
              ,(sum-ccs attm:tm-n-segments-created l2l4tm:patch-TMs) " TM segments\n"
              ,(sum-ccs attm:tm-n-synapses-created l2l4tm:patch-TMs) " TM synapses\n"
              ))))
      ))

(define exp4a '(
    [figure         .  f4a]
    [num-sequences  .   50]
    [num-objects    .    0]
    [num-features   .  100]
    [num-locations  .  100]
    [train-keys     . (TMac)]
    [test-keys      . (L4lnp L4lpa TMlnp TMlpa)]))
                                                                                            ;
(define exp5a '(
    [figure         .  f5a]
    [num-objects    .   50]
    [num-sequences  .    0]
    [num-features   .  100]
    [num-locations  .  100]
    [train-keys     . (L4ac)]
    [test-keys      . (L4lnp L4lpa TMlnp TMlpa)]))
                                                                                            ;
(define exp6 '(
    [figure          .  f6]
    [num-sequences   .  50]
    [num-objects     .  50]
    [num-features    . 500]
    [num-locations   . 100]
    [num-repetitions .  30]
    [test-keys       . (L4lpa TMlpa)]))
                                                                                            ;
(define (run figure . options)
  (let ((start (statistics)))
    (experiment
      (case figure
        [("4a")    exp4a ]
        [("5a")    exp5a ]
        [("6")     exp6  ])
      options)
    (sstats-print (sstats-difference (statistics) start))))
  
(define (option name parameter)
  ;; accept a few run options on command line
  (case name
    [("-cc")    `[column-count         . ,(string->number parameter)] ]
    [("-in")    `[intersperse-noise    . ,(string->number parameter)] ]
    [("-it")    `[interleave-training  . ,(string=? "#t"  parameter)] ]
    [("-ncc")   `[num-cortical-columns . ,(string->number parameter)] ]
    [("-nf")    `[num-features         . ,(string->number parameter)] ]
    [("-no")    `[num-objects          . ,(string->number parameter)] ]
    [("-nr")    `[num-repetitions      . ,(string->number parameter)] ]
    [("-ns")    `[num-sequences        . ,(string->number parameter)] ]
    [("-nsl")   `[no-seq-location      . ,(string=? "#t"  parameter)] ]
    [("-rsl")   `[random-seq-location  . ,(string=? "#t"  parameter)] ]
    [("-ss")    `[superimpose-sequence . ,(string=? "#t"  parameter)] ]
    [("-st")    `[settling-time        . ,(string->number parameter)] ]))
  
(define (main command-line)              ;; {String} ->
  (let process ((args (reverse command-line)) (options (list)) (parameter "#t"))
    (cond
      [(or (null? args)
           (string=? "" (car args))) ]   ;; (not run as --program)
      [(null? (cdr args))
        (run (car args) options)]        ;; last (first) element is experiment
      [(char=? #\- (string-ref (car args) 0))
        (process (cdr args) (cons (option (car args) parameter) options) "#t")]
      [else
        (process (cdr args) options (car args))])))

(random-seed! #;42 (time-second (current-time)))
        
(main (command-line-arguments))