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
  (prefix (HTM-scheme HTM-scheme algorithms l2l4tm_patch) l2l4tm:)
  (prefix (HTM-scheme HTM-scheme algorithms apical_tiebreak_pair_memory) atpm:)
  (prefix (HTM-scheme HTM-scheme algorithms apical_tiebreak_sequence_memory) atsm:)
  (prefix (HTM-scheme HTM-scheme algorithms column_pooler) cp:))
                                                                                            ;
;; CColX          = range num-cortical-columns
;; Sensation      = (FeatureSDRX, LocationSDRX)
;; Object         = vectorof Sensation
;; Sequence       = vectorof (FeatureSDRX, #f)
;; Experience     = Object|Sequence
;;   interpretation:
;;   an experience is either a succession of Sensations (feature at location)
;;   or a sequence of features (random locations will be provided)
;; ExperienceX    = range num-objects|num-sequences
;; ExperiencePool = vectorof Experience indexed by ExperienceX
;; FeaturePool    = vectorof SDR indexed by FeatureSDRX
;; LocationPool   = vectorof SDR indexed by LocationSDRX
;; SDR            = listof ColX

(define-record-type sensation (fields
  feature                                ;; SDRX
  location                               ;; SDRX
  ))
                                                                                            ;
(define-record-type om (fields           ;; Object Machine
  num-cortical-columns
  experiences                            ;; vectorof Experience
  feature-pool                           ;; vectorof SDR
  location-pool))                        ;; vectorof SDR
                                                                                            ;
(define (have-sensations-of experience   ;; Experience OM Patch Nat Number Bool [(Nat -> )]
          om patch settling-time noise learn . report)
  ;; feed each sensation of experience to the patch settling-time times, report after each sensation
  (let* ( (ns         (vector-length experience))
          (object?    (sensation-location (vector-ref experience 0)))
          (experience (if (and object? (not learn))    ;; train objects on shuffled sensations
                          (vector-sample experience ns)
                          experience))
          (ncc        (om-num-cortical-columns om)))
    (define (build-input field sx)
      (let ((fieldxs (vector-map field experience)))
        (build-vector ncc
          (if (or object? (eq? field sensation-feature))
              (lambda (ccx) 
                (+ (* ncc (vector-ref fieldxs sx)) ccx))
              (lambda _ (random (vector-length (om-location-pool om))))))))
    (do ((sx 0 (add1 sx))) ((= sx ns))
      (do ((_ 0 (add1 _))) ((= _ settling-time))
        (l2l4tm:compute patch
          (vector-refs (om-feature-pool om)  (build-input sensation-feature sx))
          (vector-refs (om-location-pool om) (build-input sensation-location sx))
          learn))
      (unless (null? report)
        ((car report) sx)))))
                                                                                            ;
(define (create-random-experiences       ;; Nat Nat Nat Nat -> ExperiencePool
          num-exps num-sensations num-features num-locations)
  ;; produce Object experiences, or Sequence experiences when num-locations is zero
  (let ((location-array (build-vector num-locations id)))
    (build-vector num-exps (lambda _
      ;; new selection of features and locations for each experience
      (let ((locations (if (positive? num-locations)
                           (vector-sample location-array num-sensations)
                           (make-vector num-sensations #f))))
        (build-vector num-sensations (lambda (sx)
            (make-sensation (random num-features)
              (vector-ref locations sx)))))))))
                                                                                            ;
(define (run-experiment args)            ;; KWargs ->
  ;; train l2l4tm-patch on object/sequence experiences, test, summarize response
  (define (get key default)
    (let ((specified (assoc key args)))
      (if specified (cdr specified) default)))
  (let* (
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
      (num-input-bits                    (get 'input-bits          20))
      (settling-time                     (get 'settling-time
                                              (if (> num-cortical-columns 1) 3 1)))
      (num-repetitions                   (get 'num-repetitions      5))
      (figure                            (get 'figure             "?"))
      (syn-perm-proximal-dec-L2          (get 'syn-perm-proximal-dec-L2 (perm 0.001)))
      (min-threshold-proximal-L2         (get 'min-threshold-proximal-L2    10))
      (sample-size-proximal-L2           (get 'sample-size-proximal-L2      15))
      (basal-predicted-segment-decrement (get 'basal-predicted-segment-decrement (perm 0.0006)))
      (num-learning-points               (get 'num-learning-points 1))
      (enable-feedback                   (get 'enable-feedback #f))
      (_                                 (random-seed! 42))
      (stats-keys                        (get 'stats-keys '(L4lpa TMlnp TMlpa)))
;; create the sequences and objects
      (generate-pattern (lambda _
                          (list-sort fx<? (vector->list 
                            (vector-sample (build-vector input-size id) num-input-bits)))))
      (generate-pool    (lambda (n) 
                          (build-vector (* n 3 #;num-cortical-columns) generate-pattern)))
      (feature-pool     (generate-pool num-features))
      (location-pool    (generate-pool num-locations))
      (sequence-experiences
        (create-random-experiences num-sequences seq-length num-features 0))
      (object-experiences
        (create-random-experiences num-objects num-points num-features num-locations))
      (sequences   
        (make-om num-cortical-columns sequence-experiences feature-pool location-pool))
      (objects                           ;; objects have same features as sequences
        (make-om num-cortical-columns object-experiences feature-pool location-pool))
;; setup experiment: multi-column multi-layer patch
      (L2-overrides  (get 'L2-overrides
                      `([syn-perm-proximal-dec       . ,syn-perm-proximal-dec-L2]
                        [min-threshold-proximal      . ,min-threshold-proximal-L2]
                        [sample-size-proximal        . ,sample-size-proximal-L2]
                        [initial-proximal-permanence . ,(perm 0.45)]
                        [syn-perm-proximal-dec       . ,(perm 0.002)])))
      (L4-overrides  (get 'L4-overrides
                      `([initial-permanence                . ,(perm 0.21)]
                        [activation-threshold              . 18]
                        [min-threshold                     . 18]
                        [basal-predicted-segment-decrement . ,basal-predicted-segment-decrement])))
      (TM-overrides  (get 'TM-overrides
                      `([basal-predicted-segment-decrement . ,basal-predicted-segment-decrement])))
      (patch 
        (l2l4tm:make-patch num-cortical-columns input-size num-input-bits
          enable-feedback L2-overrides L4-overrides TM-overrides))
      )
#;> (define (train-superimposed objects sequences . report)
      (let* (
          (sequence-order (build-list num-sequences id))
          (sequence-order (apply append (make-list num-repetitions sequence-order)))
          (ns             (length sequence-order))
          (sequence-order (vector-sample (list->vector sequence-order) ns)))
        (do ((ox 0 (add1 ox))) ((= ox num-objects))
          (let ((object-sensations (vector-ref (om-experiences objects) ox)))
            (do ((s 0 (add1 s))) ((= s num-repetitions))
              (let* (
                  (ns                  (vector-length object-sensations))
                  (object-sensations   (vector-sample object-sensations ns))
                  (sequence-id         (vector-ref sequence-order (+ (* ox num-repetitions) s)))
                  (sequence-sensations (vector-ref (om-experiences sequences) sequence-id))
                  (sxs                 (build-vector ns id)))
                (vector-for-each             ;; sensation in sensations
                  (lambda (sx)
                    (let ((obj-feature  (vector-ref (om-feature-pool objects)
                                          (sensation-feature
                                            (vector-ref object-sensations sx))))
                          (obj-location (vector-ref (om-location-pool objects)
                                          (sensation-location
                                            (vector-ref object-sensations sx))))
                          (seq-feature  (vector-ref (om-feature-pool sequences)
                                          (sensation-feature
                                            (vector-ref sequence-sensations sx)))))
                      (do ((_ 0 (add1 _))) ((= _ settling-time))
                        (l2l4tm:compute patch
                          (make-vector num-cortical-columns (union1d obj-feature seq-feature))
                          (make-vector num-cortical-columns obj-location)
                          #t))
                      (unless (null? report)
                        ((car report) sx))))
                  sxs)))))))
                                                                                            ;
#;> (define (train-all om num-repeats)   ;; OM Nat ->
      ;; train the network on all the experiences
      (vector-for-each (lambda (experience)
          (do ((_ 0 (add1 _))) ((= _ num-repeats))
            (have-sensations-of experience om patch num-learning-points 0.0 #t)
            (unless (sensation-location (vector-ref experience 0))   ;; reset TM between sequences
              (vector-for-each atsm:reset (l2l4tm:patch-TMs patch))))
          (l2l4tm:reset patch))
        (om-experiences om)))
                                                                                            ;
#;> (define (infer om ex report)         ;; OM Nat Stats
      ;; test experience, update stats after each sensation
      (have-sensations-of (vector-ref (om-experiences om) ex) 
                          om patch settling-time 0.0 #f report)
      (l2l4tm:reset patch))
                                                                                            ;
#;> (define (get-stats-for keys om ex sx)
      (define (lengths cells layer)
        (vector-map length
          (vector-map cells (layer patch))))
      (map (lambda (key)
        (case key
          [(L4lpa) (lengths atpm:get-predicted-active-cells l2l4tm:patch-L4s)]
          [(TMlnp) (lengths atsm:get-next-predicted-cells   l2l4tm:patch-TMs)]
          [(TMlpa) (lengths atpm:get-predicted-active-cells l2l4tm:patch-TMs)]
          [(L2a)   (list->vector (cp:get-active-cells (vector-ref (l2l4tm:patch-L2s patch) 0)))]      
          [(L2la)  (lengths cp:get-active-cells l2l4tm:patch-L2s)]))       
        keys))
                                                                                            ;
#;> (define (infer-all om stats)         ;; OM Stats
      (vector-for-each (lambda (ex)
        (infer om ex (lambda (sx)
          (vector-set! stats ex (append (vector-ref stats ex)
                                        (list (get-stats-for stats-keys om ex sx)))))))
        (indexes (om-experiences om))))
                                                                                            ;
;; Training:
    (case figure
      [(A6)
        (train-superimposed objects sequences)]
      [(H3b H3c)
        (train-all objects 5)]           ;; for HAC2017 figures
      [else
        (train-all objects num-repetitions)
        (train-all sequences num-repetitions)])
;; Inference:
    (let ((stats (make-vector 50 '())))
      (case figure
        [(A6)
          (let ((objectxs (vector-sample (build-vector num-objects id) 8)))
          (for-each
            (lambda (item-type ex)
              (infer
                (if (eq? item-type 'seq)  sequences  objects )
                (vector-ref objectxs ex)
                (lambda (sx)
                  (vector-set! stats 0 
                    (append (vector-ref stats 0)
                            (list (get-stats-for stats-keys #f ex sx)))))))
            '(seq obj seq obj seq seq obj seq) (build-list 8 id)))]
        [(H3b H3c)
          (infer-all objects   stats)]
        [else
          (infer-all objects   stats)
          (infer-all sequences stats)])
;; Results:
    (let ((mean (lambda (x) (if (eq? figure 'A6)  x
                   (int<- (/ x (+ num-objects num-sequences)))))))
    (case figure
      [(H3b H3c) #f]
      [else
        ;; stats :: (ExpVectorOf (SensListOf (KeyListOf (CCVectorOf Number))))
        (vector-for-each-x                   ;; for each experience
          (lambda (exp-stats ex)
            (unless (null? exp-stats)
              (for-each display (list "\n" ex ": " 
                (if (> num-cortical-columns 1)
                  (if (= (length stats-keys) 1)  exp-stats
                    (list (car exp-stats) "..." (car (reverse exp-stats))))
                  (map (lambda (sens-stats)
                          (map (lambda (v) (vector-ref v 0)) sens-stats))
                    exp-stats))))))
          stats)])
    (let ((summary-file "HTM-scheme/projects/combined_sequences/combined_sequences.data"))
      (when (file-exists? summary-file) (delete-file summary-file))
      (with-output-to-file summary-file (lambda ()
        (write [cons [list "figure" (symbol->string figure)]
          (case figure
            [(H3b H3c)
              (let ((first-object (apply append (vector-ref stats 0))))
                [cons [list "L2a"
                        (fold-left (lambda (acc v)
                            (list->vector (union1d (vector->list acc) (vector->list v))))
                          (car first-object)
                          (cdr first-object))]
                [list [list "L2ac" (vector-filter positive?
                  (fold-left (lambda (acc v)
                      (vector-for-each (lambda (x)
                          (vector-set! acc x (add1 (vector-ref acc x))))
                        v)
                      acc)
                    (make-vector (cp:number-of-cells (vector-ref (l2l4tm:patch-L2s patch) 0)) 0)
                    first-object)) ] ] ] )]
            [else
              (map                                     ;; each stats key
                (lambda (key kx)                       ;; Symbol Nat ->
                  [cons (symbol->string key)
                    [list (map mean                    ;; average
                      (vector-fold-left                ;; ..over experiences
                        (lambda (accs exp-stats)       ;; (listof Number) (SensListOf (KeyListOf (vectorof Number)))
                          (if (null? exp-stats)  accs
                            (map                       ;; ..each sensation, average columns
                              (lambda (ac sens-stats)  ;; Number (listof (vectorof Number)) -> Number
                                (+ ac (vector-average (list-ref sens-stats kx))))
                              accs exp-stats)))        ;; -> (listof Number)
                        (make-list (length (vector-ref stats 0)) 0)
                        stats)) ] ] )
                stats-keys (indexes stats-keys))])
             ] ))))))))

(define (run-experiment-4a)
  (run-experiment `(
    [figure          .  A4a]
    [num-sequences   .   50]
    [seq-length      .   10]
    [num-features    .  100]
    [trial-num       .    0]
    [num-objects     .    0]
    [num-locations   .  200]
    [num-repetitions .   10]             ;; (3 * seqLength?)
    [input-size      . 2048]
    [basal-predicted-segment-decrement . ,(perm 0.001)]
    )))
                                                                                            ;
(define (run-experiment-5a)
  (run-experiment `(
    [figure               .  A5a]
    [num-sequences        .    0]
    [seq-length           .   10]
    [num-features         .  100]
    [trial-num            .    4]
    [num-objects          .   50]
    [num-locations        .  100]
    [num-cortical-columns .    1]
    [num-learning-points  .    2]
    )))
                                                                                            ;
(define (run-experiment-6)
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
    [stats-keys      . (L4lpa TMlpa)]
    )))
                                                                                            ;
(define (run-experiment-H3b)
  (run-experiment `(
    [figure              . H3b]
    [num-cortical-columns .  1]
    [num-features         .  3]
    [num-points           . 10]
    [num-objects          . 10]
    [num-locations        . 10]
    [input-size         . 1024]
    [num-input-bits       . 20]
    [L2-overrides         . ()]
    [L4-overrides         . ()]
    [enable-feedback      . #t]
    [num-sequences        .  0]
    [num-learning-points  .  3]
    [settling-time        .  2]
    [stats-keys           . (L2a)]
    )))
                                                                                            ;
(define (run-experiment-H3c)
  (run-experiment `(
    [figure              . H3c]
    [num-cortical-columns .  3]
    [num-features         .  3]
    [num-points           . 10]
    [num-objects          . 10]
    [num-locations        . 10]
    [input-size         . 1024]
    [num-input-bits       . 20]
    [L2-overrides         . ()]
    [L4-overrides         . ()]
    [enable-feedback      . #t]
    [num-sequences        .  0]
    [num-learning-points  .  3]
    [settling-time        .  3]
    [stats-keys           . (L2a)]
    )))
                                                                                            ;