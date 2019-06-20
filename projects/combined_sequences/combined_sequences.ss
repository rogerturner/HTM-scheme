#!chezscheme

;; === HTM-scheme/projects/combined_sequences Copyright 2018 Roger Turner. ===
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

  ;; Partial replication of experiments reported in Numenta papers 
  ;; Ahmad S and Hawkins J (2017) "Untangling Sequences: Behavior vs. External Causes"
  ;; doi: 10.1101/190678, and  
  ;; Hawkins J, Ahmad S and Cui Y (2017) "A Theory of How Columns in the Neocortex 
  ;; Enable Learning the Structure of the World" doi: 10.3389/fncir.2017.00081
  ;; Translated from numenta htmresearch/.../combined_sequences.py,
  ;; combined_sequence_network_creation.py, /combined_sequence_experiment.py,
  ;; laminar_network.py, l2_l4_inference.py --
  ;; see comments in numenta code for descriptions of functions and parameters.
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.

(import
          (except (chezscheme) add1 make-list random reset)
          (HTM-scheme HTM-scheme algorithms htm_prelude)
          (HTM-scheme HTM-scheme algorithms htm_concept)
  (prefix (HTM-scheme HTM-scheme algorithms l2l4tm_patch)                    l2l4tm:)
  (prefix (HTM-scheme HTM-scheme algorithms apical_tiebreak_pair_memory)     atpm:)
  (prefix (HTM-scheme HTM-scheme algorithms apical_tiebreak_sequence_memory) atsm:)
  (prefix (HTM-scheme HTM-scheme algorithms column_pooler)                   cp:))

;; === Types ===
;; CColX          = Range num-cortical-columns
;; Sensation      = (FeatureSDRX, LocationSDRX)
;; Object         = Vectorof Sensation
;; Sequence       = Vectorof (FeatureSDRX, #f)
;; Experience     = Object|Sequence
;;   interpretation:
;;   an experience is either a succession of Sensations (feature at location)
;;   or a sequence of features (random locations will be provided)
;; ExperienceX    = Range num-objects|num-sequences
;; ExperiencePool = Vectorof Experience indexed by ExperienceX
;; FeaturePool    = Vectorof SDR indexed by FeatureSDRX
;; LocationPool   = Vectorof SDR indexed by LocationSDRX
;; SDR            = Listof ColX

(define-record-type sensation (fields
  feature                                ;; SDRX
  location                               ;; SDRX
  ))
                                                                                            ;
(define (run-experiment exp-args run-args) ;; KWargs KWargs -> 
  ;; train l2l4tm-patch on object/sequence experiences, test, summarize response
  (let* (
      (args (append run-args exp-args))
      (get  (lambda (key default)
              (let ((specified (assoc key args)))
                (if specified (cdr specified) default))))
;; experiment parameters
      (num-cortical-columns              (get 'num-cortical-columns 1))
      (num-objects                       (get 'num-objects         10))
      (num-sequences                     (get 'num-sequences       10))
      (num-features                      (get 'num-features        10))
      (seq-length                        (get 'seq-length          10))
      (num-points                        (get 'num-points          10))
      (trial-num                         (get 'trial-num           42))
      (input-size                        (get 'input-size        1024))
      (num-locations                     (get 'num-locations   100000))
      (num-input-bits                    (get 'num-input-bits      20))
      (settling-time                     (get 'settling-time        1))
      (num-repetitions                   (get 'num-repetitions      5))
      (figure                            (get 'figure              'z))
      (syn-perm-proximal-dec-L2          (get 'syn-perm-proximal-dec-L2 (perm 0.001)))
      (min-threshold-proximal-L2         (get 'min-threshold-proximal-L2 10))
      (sample-size-proximal-L2           (get 'sample-size-proximal-L2   15))
      (sample-size-distal-L2             (get 'sample-size-distal-L2     20))
      (basal-predicted-segment-decrement (get 'basal-predicted-segment-decrement (perm 0.001)))
      (num-learning-points               (get 'num-learning-points  1))
      (num-cells-per-column              (get 'num-cells-per-column 16))
      (enable-feedback                   (get 'enable-feedback     #f))
      (online-learning                   (get 'online-learning     #f))
      (random-sequence-location          (get 'random-sequence-location  #t))
      (display-timing                    (get 'display-timing      #t))
      (train-keys                        (get 'train-keys '()))
      (test-keys                         (get 'test-keys  '()))
      (train-keys (case figure
                    [(H3b H3c H4b) (take num-cortical-columns train-keys)]
                    [else train-keys]))
      (test-keys  (case figure
                    [(H3b H3c H4b) (take num-cortical-columns test-keys)]
                    [else test-keys]))
      #;(_          (random-seed! 42))
;; create the sequences and objects
#;>   (random-sdr   (lambda (size)       ;; Nat -> (X -> SDR)
                    ;; produce function to generate SDR of size
                      (lambda _
                        (list-sort fx<? (vector->list 
                          (vector-sample (build-vector size id) num-input-bits))))))
      (feature-pool
        (build-vector (* num-features  num-cortical-columns) (random-sdr input-size)))
      (location-pool
        (build-vector (* num-locations num-cortical-columns) (random-sdr (* input-size num-cells-per-column))))
#;>   (create-random-experiences         ;; Nat Nat Nat -> Experiences
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
      (sequences (create-random-experiences num-sequences seq-length 0))
      (objects   (create-random-experiences num-objects num-points num-locations))
;; setup experiment
      (L2-overrides  (get 'L2-overrides
                      `([syn-perm-proximal-dec       . ,syn-perm-proximal-dec-L2]
                        [min-threshold-proximal      . ,min-threshold-proximal-L2]
                        [sample-size-proximal        . ,sample-size-proximal-L2]
                        [sample-size-distal          . ,sample-size-distal-L2]
                        [initial-proximal-permanence . ,(perm 0.45)]
                        [syn-perm-proximal-dec       . ,(perm 0.002)]
                        [predicted-inhibition-threshold . ,(int<- (* 0.8 num-input-bits))]
                        [online-learning             . ,online-learning])))
      (TM-overrides  (get 'TM-overrides
                      `([basal-predicted-segment-decrement . ,basal-predicted-segment-decrement])))
      (L4-overrides  (get 'L4-overrides
                      `([initial-permanence                . ,(perm 0.21)]
                        [activation-threshold              . 18]
                        [min-threshold                     . 18]
                        [basal-predicted-segment-decrement . ,basal-predicted-segment-decrement])))
      (patch 
        (l2l4tm:make-patch num-cortical-columns input-size num-input-bits
          enable-feedback L2-overrides L4-overrides TM-overrides))
      (start-time     (cpu-time))
      )
#;> (define (train-superimposed)         ;; -> 
      ;; train objects with feature superimposed with feature from random sequences
      (let* ( (sequence-order (apply append (make-list num-repetitions 
                                              (build-list num-sequences id))))
              (ns             (length sequence-order))
              (sequence-order (vector-sample (list->vector sequence-order) ns)))
        (do ((ox 0 (add1 ox))) ((= ox num-objects))
          (let ((object-sensations (vector-ref objects ox)))
            (do ((s 0 (add1 s))) ((= s num-repetitions))
              (let* ( (ns          (vector-length object-sensations))
                      (object-ss   (vector-sample object-sensations ns))
                      (sequence-id (vector-ref sequence-order (+ (* ox num-repetitions) s)))
                      (sequence-ss (vector-ref sequences sequence-id)))
                (vector-for-each         ;; sensation in sensations
                  (lambda (sx)
                    (let ((obj-feature (union1d
                            (vector-ref feature-pool
                              (sensation-feature (vector-ref object-ss sx)))
                            (vector-ref feature-pool
                              (sensation-feature (vector-ref sequence-ss sx)))))
                          (obj-location
                            (vector-ref location-pool
                              (sensation-location (vector-ref object-ss sx)))))
                      (do ((_ 0 (add1 _))) ((= _ settling-time))
                        (l2l4tm:compute patch
                          (make-vector num-cortical-columns obj-feature)
                          (make-vector num-cortical-columns obj-location)
                          #t))))
                  (build-vector ns id))))))))
                                                                                            ;
#;> (define (infer-switching stats)      ;; Stats -> 
      ;; test objects & sequences after superimposed training
      (let ((objectxs (vector-sample (build-vector num-objects id) 8)))
        (for-each
          (lambda (item-type ex)
            (infer
              (if (eq? item-type 'seq)  sequences  objects )
              (vector-ref objectxs ex)
              (lambda (sx)
                (vector-set! stats 0 
                  (append (vector-ref stats 0)
                          (list (get-stats-for test-keys)))))))
          '(seq obj seq obj seq seq obj seq) (build-list 8 id))))
                                                                                            ;
#;> (define (have-sensations-of          ;; Experience Bool [(Nat -> )] ->
              experience learn . report)
      ;; feed each sensation of experience settling-time times, report after each sensation
      ;; (replaces both learnObjects and infer)
      (let* ( (object?    (sensation-location (vector-ref experience 0)))
              (ns         (vector-length experience))
              (experience (if (and object? (not learn))    ;; infer objects on shuffled sensations
                              (vector-sample experience ns)
                              experience)))
        (define (build-input field sx)
          (let ((fieldxs (vector-map field experience)))
            (build-vector num-cortical-columns
              (if (or object? (eq? field sensation-feature))
                  (lambda (ccx) 
                    (+ (* num-cortical-columns (vector-ref fieldxs sx)) ccx))
                  (lambda _ 
                    (if random-sequence-location
                      (random (vector-length location-pool))
                      sx) )))))
        (do ((sx 0 (add1 sx))) ((= sx ns))
          (do ((_ 0 (add1 _))) ((= _ settling-time))
            (l2l4tm:compute patch
              (vector-refs feature-pool  (build-input sensation-feature sx))
              (vector-refs location-pool (build-input sensation-location sx))
              (or learn online-learning)))
          (unless (null? report)
            ((car report) sx)))))
                                                                                            ;
#;> (define (train-all es stats)  ;; Experiences Nat Stats -> 
      ;; train the network on all the experiences (replaces trainSequences and trainObjects)
      (vector-for-each 
        (lambda (e ex)
          (do ((_ 0 (add1 _))) ((= _ num-repetitions))
            (have-sensations-of e #t)
            (unless (sensation-location (vector-ref e 0))
              ;; reset TM between presentations of sequence (num-learning-points = 1 for sequences)
              (vector-for-each atsm:reset (l2l4tm:patch-TMs patch))))
          (vector-set! stats ex 
            (append (vector-ref stats ex)
                    (list (get-stats-for train-keys))))
          (l2l4tm:reset patch))
        es (indexes es)))
                                                                                            ;
#;> (define (infer es ex report)         ;; Experiences Nat (Nat -> )
      ;; test experience, update stats after each sensation
      (have-sensations-of (vector-ref es ex) #f report)
      (l2l4tm:reset patch))              ;; reset defaults to true for infer
                                                                                            ;
#;> (define (get-stats-for keys)         ;; {Symbol} -> {Nat}
      (define (lengths cells layer)
        (vector-map length
          (vector-map cells (layer patch))))
      (map (lambda (key)
        (case key
          [(L4lpa) (lengths atpm:get-predicted-active-cells l2l4tm:patch-L4s)]
          [(TMlnp) (lengths atsm:get-next-predicted-cells   l2l4tm:patch-TMs)]
          [(TMlpa) (lengths atsm:get-predicted-active-cells l2l4tm:patch-TMs)]
          [(L2r)   (cp:get-active-cells (vector-ref (l2l4tm:patch-L2s patch) 0))]      
          [(L2r1)  (cp:get-active-cells (vector-ref (l2l4tm:patch-L2s patch) 1))]      
          [(L2r2)  (cp:get-active-cells (vector-ref (l2l4tm:patch-L2s patch) 2))]      
          [(L2a)   (cp:get-active-cells (vector-ref (l2l4tm:patch-L2s patch) 0))]      
          [(L2a1)  (cp:get-active-cells (vector-ref (l2l4tm:patch-L2s patch) 1))]      
          [(L2a2)  (cp:get-active-cells (vector-ref (l2l4tm:patch-L2s patch) 2))]))
        keys))
                                                                                            ;
#;> (define (infer-all es stats)         ;; Experiences Stats
      ;; test each experience and save stats for each sensation
      (vector-for-each (lambda (ex)
          (infer es ex 
            (lambda (sx)
              (vector-set! stats ex 
                (append (vector-ref stats ex)
                        (list (get-stats-for test-keys)))))))
        (indexes es)))
                                                                                            ;
    (let ((stats (make-vector (max num-objects num-sequences) '())))
;; train the network
      (case figure
        [(A6)
          (train-superimposed)]
        [(H3b H3c H4b)
          (train-all objects stats)]
        [else
          (train-all objects   stats)
          (train-all sequences stats)])
;; run inference
      (case figure
        [(A6)
          (infer-switching stats) ]
        [(H3b H3c H4b)
          (infer objects 0 (lambda (sx)  ;; just first object
              (vector-set! stats 0 
                (append (vector-ref stats 0)
                        (list (get-stats-for test-keys)))))) ]
        [else
          (infer-all objects   stats)
          (infer-all sequences stats) ] )
;; display timing
      (when display-timing
        (let ((tenths (fxdiv (+ (- (cpu-time) start-time) 50) 100)))
          (for-each display `(
              ,(symbol->string figure) " " ,@run-args #\newline
              ,(fxdiv tenths 10) "." ,(fxmod tenths 10) "s cpu time\n"
              ,(atpm:get-n-segments-created (vector-ref (l2l4tm:patch-L4s patch) 0)) " L4 segments\n"
              ,(atpm:get-n-synapses-created (vector-ref (l2l4tm:patch-L4s patch) 0)) " L4 synapses\n"
              ,(atsm:get-n-segments-created (vector-ref (l2l4tm:patch-TMs patch) 0)) " TM segments\n"
              ,(atsm:get-n-synapses-created (vector-ref (l2l4tm:patch-TMs patch) 0)) " TM synapses\n"
              ))))
;; compute statistics needed for plots
    (let ((mean (lambda (x) (if (eq? figure 'A6)  x
                   (int<- (/ x (+ num-objects num-sequences)))))))
    (case figure
      [(H3b H3c H4b) #f]
      [else
        ;; stats :: (ExpVectorOf (SensListOf (KeyListOf (CCVectorOf Number))))
        #f #;(vector-for-each-x               ;; for each experience (diagnostic output)
          (lambda (exp-stats ex)
            (unless (null? exp-stats)
              (for-each display (list "\n" ex ": " 
                (if (> num-cortical-columns 1)
                  (if (= (length test-keys) 1)  exp-stats
                    (list (car exp-stats) "..." (car (reverse exp-stats))))
                  (map (lambda (sens-stats)
                          (map (lambda (v) (vector-ref v 0)) sens-stats))
                    exp-stats))))))
          stats)])
    (let ((summary-file "HTM-scheme/projects/combined_sequences/combined_sequences.data"))
      (with-output-to-file summary-file (lambda ()
        (write [cons* 
          [list "figure" (symbol->string figure)]
          [list "using"  run-args]
          (case figure
            [(H3b H3c H4b)
              (let ((first-object (vector-ref stats 0)))  ;; (SensListOf (KeyListOf (CCVectorOf Number)))
                (append
                  (map                   ;; cell indexes for each stat
                    (lambda (key kx)
                      [list (symbol->string key)
                            (list-ref (car first-object) kx)])
                    train-keys (indexes train-keys))
                  (map                   ;; cell indexes for each stat
                    (lambda (key kx)
                      [list (symbol->string key)
                                         ;; fold over sensations
                            (fold-left (lambda (acc v)
                                         (union1d acc (list-ref v kx)))
                              '()
                              first-object)])
                    test-keys (indexes test-keys))
                  (map                   ;; number of sensations/cell for each stat
                    (lambda (key kx)
                      [list (string-append (symbol->string key) "c")
                            (filter positive? (vector->list
                                (fold-left (lambda (acc v)
                                    (for-each (lambda (x)
                                        (vector-set! acc x (add1 (vector-ref acc x))))
                                      (list-ref v kx))
                                    acc)
                                  (make-vector (cp:number-of-cells (vector-ref (l2l4tm:patch-L2s patch) 0)) 0)
                                  first-object))) ])
                    test-keys (indexes test-keys))))]
            [else                        ;; AH2017 figures
              (map                                     ;; each stats key
                (lambda (key kx)                       ;; Symbol Nat ->
                  [cons (symbol->string key)
                    [list (map mean                    ;; average of objects/sequences for 4a/5a
                      (vector-fold-left                ;; ..over experiences
                        (lambda (accs exp-stats)       ;; (listof Number) (SensListOf (KeyListOf (vectorof Number)))
                          (if (null? exp-stats)  accs
                            (map                       ;; ..each sensation, average columns
                              (lambda (ac sens-stats)  ;; Number (listof (vectorof Number)) -> Number
                                (if (null? sens-stats)  ac
                                  (+ ac (vector-average (list-ref sens-stats kx)))))
                              accs exp-stats)))        ;; -> (listof Number)
                        (make-list (length (vector-ref stats 0)) 0)
                        stats)) ] ] )
                test-keys (indexes test-keys))])
             ] )) 'truncate ))))))

(define (run-experiment-A4a . args)
  (run-experiment `(
    [figure          .  A4a]
    [num-sequences   .   50]
    [seq-length      .   10]
    [num-features    .  100]
    [trial-num       .    0]
    [num-objects     .    0]
    [num-locations   .  200]
    [num-repetitions .   30]
    [input-size      . 2048]
    [basal-predicted-segment-decrement . ,(perm 0.001)]
    [test-keys       . (L4lpa TMlnp TMlpa)])
    args))
                                                                                            ;
(define (run-experiment-A5a . args)
  (run-experiment `(
    [figure               . A5a]
    [num-sequences        .   0]
    [seq-length           .  10]
    [num-features         . 100]
    [trial-num            .   4]
    [num-objects          .  50]
    [num-locations        . 100]
    [test-keys            . (L4lpa TMlnp TMlpa)])
    args))
                                                                                            ;
(define (run-experiment-A6  . args)
  (run-experiment `(
    [figure          .  A6]
    [num-sequences   .  50]
    [seq-length      .  10]
    [num-objects     .  50]
    [num-features    . 500]
    [trial-num       .   8]
    [num-locations   . 100]
    [settling-time   .   1]
    [num-repetitions .  30]
    [basal-predicted-segment-decrement . ,(perm 0.001)]
    [test-keys       . (L4lpa TMlpa)])
    args))
                                                                                            ;
(define (run-experiment-H3b . args)
  (run-experiment `(
    [figure              . H3b]
    [input-size         . 1024]          ;; l2_pooling/convergence_activity.py
    [num-input-bits       . 20]
    [num-learning-points  .  4]          ;; ?
    [num-cortical-columns .  1]
    [num-features         .  3]
    [num-points           . 10]
    [num-objects          . 10]
    [num-locations        . 10]
    [L4-overrides         . ()]
    [enable-feedback      . #t]
    [num-sequences        .  0]
    [settling-time        .  1]
    [train-keys           . (L2r)]
    [test-keys            . (L2a)])
    args))
  ;; Other settings (cf Hawkins Ahmed Cui 2017):
  ;;   (run-experiment-H3b '[input-size . 150] '[num-input-bits . 10])
  ;;   (run-experiment-H3b '[enable-feedback . #f] '[num-learning-points . 5])
  ;;   (run-experiment-H3b '[num-features . 30] '[num-locations . 30] '[num-objects . 30])
                                                                                            ;
(define (run-experiment-H3c . args)
  (run-experiment `(
    [figure              . H3c]
    [input-size         . 1024]          ;; l2_pooling/convergence_activity.py
    [num-input-bits       . 20]
    [num-learning-points  .  4]          ;; ?
    [num-cortical-columns .  3]
    [num-features         .  3]
    [num-points           . 10]
    [num-objects          . 10]
    [num-locations        . 10]
    [L4-overrides         . ()]
    [enable-feedback      . #t]
    [num-sequences        .  0]
    [settling-time        .  2]
    [train-keys           . (L2r L2r1 L2r2)]
    [test-keys            . (L2a L2a1 L2a2)])
    args))
  ;; Other settings (cf Hawkins Ahmed Cui 2017):
  ;;   (run-experiment-H3c '[input-size . 150] '[num-input-bits . 10])

(define (run-experiment-H4b . args)
  (run-experiment `(
    [figure              . H4b]
    [input-size          . 150]
    [num-input-bits       . 10]
    [num-learning-points  .  4]
    [num-cortical-columns .  1]
    [num-features         .  5]
    [num-points           . 10]
    [num-objects          . 100]
    [num-locations        . 100]
    [L4-overrides         . ()]
    [enable-feedback      . #t]
    [num-sequences        .  0]
    [settling-time        .  2]
    [train-keys           . (L2r L2r1 L2r2)]
    [test-keys            . (L2a L2a1 L2a2)])
    args))

(for-each
  (lambda (expt)
    (case expt
      [("A4a")   (run-experiment-A4a) ]
      [("A4a10") (run-experiment-A4a '[num-sequences . 10]) ]
      [("A5a")   (run-experiment-A5a) ]
      [("A6")    (run-experiment-A6)  ]
      [("A610")  (run-experiment-A6  '[num-objects . 10] '[num-sequences . 10]) ]
      [("H3b")   (run-experiment-H3b '[input-size . 150] '[num-input-bits . 10]) ]
      [("H3c")   (run-experiment-H3c '[input-size . 150] '[num-input-bits . 10]) ]
      [("H4b")   (run-experiment-H4b) ]
      )
    (display-statistics))
  (command-line-arguments))

