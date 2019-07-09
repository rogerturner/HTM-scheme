#!chezscheme

;; === HTM-scheme/projects/untangling_sequences Copyright 2019 Roger Turner. ===
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
  ;; Ahmad S and Hawkins J (2017) "Untangling Sequences: Behavior vs. External Causes"
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
  (let* (
      (args (append run-args exp-args))
      (get  (lambda (key default)
              (let ((specified (assoc key args)))
                (if specified (cdr specified) default))))
;; experiment parameters
      (num-cortical-columns              (get 'num-cortical-columns  1))
      (column-count                      (get 'column-count        512))
      (num-input-bits                    (get 'num-input-bits       20))
      (activation-threshold              (get 'activation-threshold 18))
      (min-threshold                     (get 'min-threshold        18))
      (num-objects                       (get 'num-objects          50))
      (num-sequences                     (get 'num-sequences        50))
      (seq-length                        (get 'seq-length           10))
      (num-points                        (get 'num-points           10))
      (num-features                      (get 'num-features        100))
      (num-locations                     (get 'num-locations       100))
      (settling-time                     (get 'settling-time         1))
      (num-repetitions                   (get 'num-repetitions       5))
      (figure                            (get 'figure               'z))
      (basal-predicted-segment-decrement (get 'basal-predicted-segment-decrement (perm 0.001)))
      (num-cells-per-column              (get 'num-cells-per-column 10))
      (enable-feedback                   (get 'enable-feedback      #t))
      (online-learning                   (get 'online-learning      #f))
      (random-sequence-location          (get 'random-sequence-location  #t))
      (intersperse-noise                 (get 'intersperse-noise     0))
      (display-timing                    (get 'display-timing       #t))
      (train-keys                        (get 'train-keys          '()))
      (test-keys                         (get 'test-keys           '()))
;; create the sequences and objects
      (random-sdr 
        (lambda (size bits)         ;; Nat -> (X -> SDR)
          ;; produce function to generate SDR of size with sparsity bits/size
          (let ((range (build-vector size id)))
            (lambda _ (list-sort fx<? (vector->list 
                (vector-sample range bits)))))))
      (feature-pool                      ;; Vectorof SDR indexed by FeatureSDRX
        (build-vector num-features  (random-sdr column-count num-input-bits)))
      (location-pool                     ;; Vectorof SDR indexed by LocationSDRX
        (build-vector num-locations (random-sdr (* num-cells-per-column column-count) (* 2 num-input-bits))))
      (random-locations                  ;; Vectorof SDR indexed by LocationSDRX
        (build-vector 1000 (random-sdr (* num-cells-per-column column-count) (* 2 num-input-bits))))
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
;; setup network
      (L2-overrides  (get 'L2-overrides
                      `([cell-count . ,(int<- (/ (* column-count num-cells-per-column) 2))]
                        [sample-size-proximal           . ,(int<- (* 0.8 num-input-bits))]
                        [min-threshold-proximal         . ,(int<- (* 0.5 num-input-bits))]
                        [predicted-inhibition-threshold . ,(int<- (* 0.8 num-input-bits))]
                        [sample-size-distal             . ,num-input-bits]
                        [activation-threshold-distal    . ,(int<- (* 0.6 num-input-bits))]
                        [sdr-size                       . ,num-input-bits]
                        [syn-perm-proximal-dec          . ,(perm 0.002)]
                        [online-learning                . ,online-learning])))
      (TM-overrides  (get 'TM-overrides
                      `([initial-permanence                 . ,(perm 0.41)]
                        [activation-threshold               . ,activation-threshold]
                        [min-threshold                      . ,min-threshold]
                        [permanence-decrement               . ,(perm 0.02)]
                        [apical-predicted-segment-decrement . ,(perm 0.001)]
                        [basal-predicted-segment-decrement  . ,basal-predicted-segment-decrement])))
      (L4-overrides  (get 'L4-overrides
                      `([initial-permanence                 . ,(perm 0.41)]
                        [activation-threshold               . ,activation-threshold]
                        [min-threshold                      . ,min-threshold]
                        [permanence-decrement               . ,(perm 0.02)]
                        [apical-predicted-segment-decrement . ,(perm 0.001)]
                        [basal-predicted-segment-decrement  . ,basal-predicted-segment-decrement])))
      (patch 
        (l2l4tm:make-patch num-cortical-columns column-count num-input-bits
          enable-feedback L2-overrides L4-overrides TM-overrides))
                      ;; (ExpVectorOf (SensVectorOf (KeyListOf (CCVectorOf Number|{Nat}))))
      (stats          (build-vector (max num-objects num-sequences)
                          (lambda _ 
                            (build-vector (if (eq? figure 'f6)  (* 8 num-points)
                                              (max num-points seq-length))
                                        (lambda _ (list))))))
      (start-time     (cpu-time))
      )

#;> (define (have-sensations-of          ;; Experience Bool [(Nat -> )] ->
              experience learn mix . report)
      ;; feed each sensation of experience settling-time times, report after each sensation
      (let* ( 
          (object?    (sensation-location (vector-ref experience 0)))
          (ns         (vector-length experience))
          (sensations (if (and object? (not learn))  ;; infer objects on shuffled sensations
                          (vector-sample experience ns)
                          experience)))
        (do ((sx 0 (add1 sx))) ((= sx ns))
          (let* ( 
              ;; for each sensation of object, ccs get different selection of feature+location
              ;; for sequences, all ccs get same feature
              (select (if learn
                        (make-vector num-cortical-columns sx)
                        (vector-sample (indexes sensations) num-cortical-columns)))
              (features  (vector-refs feature-pool
                (if object?
                  (build-vector num-cortical-columns (lambda (ccx)
                    (sensation-feature (vector-ref sensations (vector-ref select ccx)))))
                  (build-vector num-cortical-columns (lambda _
                    (sensation-feature (vector-ref sensations sx)))))))
              #;(features
                (build-vector num-cortical-columns (lambda (ccx)
                  (if object?
                    (mix (vector-ref feature-pool
                           (sensation-feature (vector-ref sensations (vector-ref select ccx))))
                         sx)
                    (vector-ref feature-pool (sensation-feature (vector-ref sensations sx)))))))
              (locations
                (if object?
                  (vector-refs location-pool
                    (build-vector num-cortical-columns (lambda (ccx)
                      (sensation-location (vector-ref sensations (vector-ref select ccx))))))
                  (build-vector num-cortical-columns (lambda _ 
                      (vector-ref random-locations (random 1000))))
                  #;(build-vector num-cortical-columns (lambda _
                      (if random-sequence-location
                        (random (vector-length location-pool))
                        sx))))))
            (do ((_ 0 (add1 _))) ((= _ settling-time))
              (l2l4tm:compute patch features locations (or learn online-learning)))
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
#;> (define (train-all es mix)           ;; Experiences (Feature -> Feature) ->
      ;; train on each e num-repetitions times
      (vector-for-each 
        (lambda (e ex)
          (do ((_ 0 (add1 _))) ((>= _ (- num-repetitions 1)))
            (have-sensations-of e #t mix)
            (do ((_ 0 (add1 _))) ((= _ intersperse-noise))
              (have-sensations-of (vector-ref (create-random-experiences 1 seq-length 0) 0) #t mix))
            (unless (sensation-location (vector-ref e 0))
              ;; reset TM between presentations of sequence
              (vector-for-each attm:reset (l2l4tm:patch-TMs patch))))
          (have-sensations-of e #t mix   ;; last training repetition: save stats
            (lambda (sx features)
              (vector-set! (vector-ref stats ex) sx 
                (append (get-stats-for train-keys) (list (cons 'feat features)) '()))))
          (l2l4tm:reset patch))
        es (indexes es)))
                                                                                            ;
#;> (define (train-all-interleaved)      ;; ->
      ;; train alternating objects & sequences, repeat num-repetitions times
      (do ((_ 0 (add1 _))) ((= _ num-repetitions))
        (let ((start-time (cpu-time)))
          (vector-for-each 
            (lambda (o s)
                (have-sensations-of o #t unmixed)
                (have-sensations-of s #t unmixed)
                (vector-for-each attm:reset (l2l4tm:patch-TMs patch)))
            objects sequences)
          (let ((secs (fxdiv (+ (- (cpu-time) start-time) 500) 1000)))
            (display secs) (display " ") (flush-output-port (current-output-port))))))
                                                                                            ;
#;> (define (infer-all es mix)           ;; Experiences (Feature -> Feature) ->
      ;; test each experience and save stats for each sensation
      (vector-for-each (lambda (ex)
          (have-sensations-of (vector-ref es ex) #f mix
            (lambda (sx features)
              (vector-set! (vector-ref stats ex) sx
                (append (get-stats-for test-keys)
                        (vector-ref (vector-ref stats ex) sx)
                        '()))))
          (l2l4tm:reset patch))
        (indexes es)))
                                                                                            ;
#;> (define (infer-switching)            ;; ->
      ;; test mix of objects & sequences
      (let ((objectxs (vector-sample (build-vector num-objects id) 8)))
        (for-each
          (lambda (item-type ix)
            (have-sensations-of
              (vector-ref
                (if (eq? item-type 'seq) sequences  objects) (vector-ref objectxs ix))
              #f unmixed
              (lambda (sx features)
                (vector-set! (vector-ref stats 0) (+ (* 10 ix) sx)
                  (get-stats-for test-keys)))))
          '(seq obj seq obj seq seq obj seq) (build-list 8 id))))
                                                                                            ;
#;> (define (unmixed feature sx)
      feature)
                                                                                            ;
#;> (define (mix object-feature sx)
      (let ((sensations (vector-ref sequences (random num-sequences))))
        (unique! fx=? 
          (union1d object-feature
                   (vector-ref feature-pool
                     (sensation-feature (vector-ref sensations sx)))))))

;; run experiment: 1 train network
    (case figure
      [(f6)
        (train-all-interleaved)]
      [else
        (train-all objects   unmixed)
        (train-all sequences unmixed)
        ])
;;                 2 reset and run inference
    (l2l4tm:reset patch)
    (case figure
      [(f6)
        (infer-switching) ]
      [else
        (infer-all objects   unmixed)
        (infer-all sequences unmixed) ] )
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
                              [(TMpa)  (add-overlap-with 'TMac)]  #;(vector (attm:cols-from-cells
                                                (vector-ref (l2l4tm:patch-TMs patch) 0)
                                                (vector-ref (cdr stats-for-key) 0)))
                              [(TMlnp TMlpa L4lnp L4lpa)
                                (+ ac #;(apply max (vector->list (cdr stats-for-key))) 
                                  (vector-average (cdr stats-for-key)))]
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

(define AH2017 '(
    [column-count         . 512]
    [num-input-bits       .  20]
    [activation-threshold .  18]
    [min-threshold        .  18]))
                                                                                            ;
(define MC150 '(
    [column-count         . 150]
    [num-input-bits       .   6]
    [activation-threshold .   4]
    [min-threshold        .   4]))
                                                                                            ;
(define exp4a '(
    [figure         .  f4a]
    [num-sequences  .    5]
    [num-objects    .    0]
    [num-features   .   10]
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
    [num-features    .  50]
    [num-locations   .  50]
    [num-repetitions .  30]
    [test-keys       . (L4lpa TMlpa)]))
                                                                                            ;
  (if (pair? (command-line-arguments))   
    (random-seed! (time-second (current-time)))
    (random-seed! 42))

(define (run exp config . args)
  (experiment (append exp config) args))

(for-each (lambda (expt)                 ;; when run as program, not in repl
    (case expt
      [("4a")    (run exp4a AH2017) ]
      [("5a")    (run exp5a AH2017) ]
      [("6")     (run exp6  AH2017) ]
      [("4am")   (run exp4a MC150 ) ]
      [("4am3")  (run exp4a MC150  '[num-cortical-columns . 3] ) ]
      [("4am7")  (run exp4a MC150  '[num-cortical-columns . 7] ) ]
      [("5am")   (run exp5a MC150 ) ]
      [("5am3")  (run exp5a MC150  '[num-cortical-columns . 3] ) ]
      [("5am7")  (run exp5a MC150  '[num-cortical-columns . 7] ) ]
      [("6m")    (run exp6  MC150 ) ]
      [("6m3")   (run exp6  MC150  '[num-cortical-columns . 3] ) ]
      [("6m7")   (run exp6  MC150  '[num-cortical-columns . 7] ) ]
      )
    (display-statistics))
  (command-line-arguments))
