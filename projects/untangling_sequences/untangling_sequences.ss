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
  ;;
  ;; "Remember that all models are wrong; the practical question is how wrong
  ;;  do they have to be to not be useful" [George Box]

(import
          (except (chezscheme) add1 make-list random reset)
          (HTM-scheme HTM-scheme algorithms htm_prelude)
          (HTM-scheme HTM-scheme algorithms htm_concept)
  (prefix (HTM-scheme HTM-scheme algorithms l2_l4_patch)  l2l4:))

;; === Types (see htm_concept.ss) ===
;; SDR            = Listof Nat
;; CColX          = Range num-cortical-columns
;; FeatureSDRX    = Range num-features
;; LocationSDRX   = Range num-locations
;; Sensation      = (FeatureSDRX, LocationSDRX|#f)
;; Object         = Vectorof (FeatureSDRX, LocationSDRX)
;; Sequence       = Vectorof (FeatureSDRX, #f)
;; Experience     = Object|Sequence

(define-record-type sensation (fields
  feature                                ;; SDRX
  location))                             ;; SDRX|#f

(define (experiment exp-args run-args)   ;; KWargs KWargs ->
  ;; train l2l4 patch on object/sequence experiences, test, summarize response
  (let* [
    (args (append run-args exp-args))    ;; run-args override exp-args
    (get  (lambda (key default)
            (let ((specified (assoc key args)))
              (if specified (cdr specified) default))))
                                                                                            ;
;; model and learning parameter defaults (overridable by experiment/run parameters)
    (num-cortical-columns         (get 'num-cortical-columns          1))
    (num-minicolumns              (get 'column-count                150))
    (num-cells-per-l4pop          (get 'num-cells-per-l4pop          10))
    (num-input-bits               (get 'num-input-bits (min 20 (max 4 
                                          (int<- (* num-minicolumns 0.03))))))
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
    (online-learning              (get 'online-learning      #f))
                                                                                            ;
;; setup network: override algorithm default parameters (scaled to values above)
    (l2-sample-size-proximal      (* 2 num-input-bits))
    (l2-sample-size-distal        (* 5 num-input-bits))
    (l2-overrides  (get 'l2-overrides `(
      [initial-proximal-permanence    . ,initial-permanence]
      [connected-permanence-proximal  . ,connected-permanence]
      [initial-distal-permanence      . ,initial-permanence]
      [connected-permanence-distal    . ,connected-permanence]
      [sample-size-proximal           . ,l2-sample-size-proximal]
      [min-threshold-proximal         . ,(max 4 (int<- (* 0.8 l2-sample-size-proximal)))]
      [predicted-inhibition-threshold . ,l2-sample-size-proximal]
      [sample-size-distal             . ,l2-sample-size-distal]
      [activation-threshold-distal    . ,(int<- (* 0.8 l2-sample-size-distal))]
      [sdr-size                       . ,(* 5 num-input-bits)]
      [online-learning                . ,online-learning])))
    (l4-overrides  (get 'l4-overrides `(
      [initial-permanence                 . ,initial-permanence]
      [connected-permanence               . ,connected-permanence]
      [permanence-increment               . ,permanence-increment]
      [permanence-decrement               . ,permanence-decrement]
      [activation-threshold               . ,activation-threshold]
      [min-threshold                      . ,min-threshold]
      [reduced-basal-threshold            . ,reduced-basal-threshold]
      [apical-predicted-segment-decrement . ,predicted-segment-decrement]
      [basal-predicted-segment-decrement  . ,predicted-segment-decrement])))
    (external-input-size (int<- (* 1.5 num-minicolumns num-cells-per-l4pop)))
    (patch 
      (l2l4:make-patch num-cortical-columns num-minicolumns num-cells-per-l4pop num-input-bits
        external-input-size l2-overrides l4-overrides))
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
    (random-seq-location               (get 'random-seq-location  #f))
    (location-per-sequence             (get 'location-per-sequence #f))
    (intersperse-noise                 (get 'intersperse-noise     0))
    (superimpose-sequence              (get 'superimpose-sequence #f))
    (interleave-training               (get 'interleave-training  #f))
                                                                                            ;
;; create the sequences and objects
    (random-sdr (lambda (size bits)    ;; Nat -> (X -> SDR)
      ;; produce function to generate SDR of size with sparsity bits/size
      (let ((range (build-vector size id)))
        (lambda () (list-sort fx<? (vector->list 
            (vector-sample range bits)))))))
    (build-distinct (lambda (n ->sdr)  ;; Nat (-> SDR) -> (vectorof SDR)
      ;; produce vector of SDRs with minimum overlaps (no more than 3 bits)
      (let ((v (build-vector n (lambda _ (->sdr)))))
        (do ((i 1 (add1 i))) ((fx=? i n))
          (call/cc (lambda (break)
              (let try-again ((tries 0) (overlap 0))
                (when (fx>? overlap 3) (error #f "can't generate distinct features or locations"))
                (if (fx>? tries 100)
                  (try-again 0 (add1 overlap))
                  (let ((vi (vector-ref v i)))
                    (do ((j (fx- i 1) (fx- j 1))) ((fxnegative? j))
                      (when (fx>? (length (intersect1d (vector-ref v j) vi)) overlap)
                        (vector-set! v i (->sdr))
                        (break (try-again (add1 tries) overlap))))))))))
        v)))
    (feature-pool                      ;; Vectorof SDR indexed by FeatureSDRX
      (build-distinct num-features  (random-sdr num-minicolumns num-input-bits)))
    (location-pool                     ;; Vectorof SDR
      (build-distinct (* 2 num-locations) 
                      (random-sdr external-input-size (int<- (* 2.5 num-input-bits)))))
    (feature-locations                 ;; Vectorof SDR indexed by LocationSDRX
      (build-vector num-locations (lambda (i) (vector-ref location-pool i))))
    (random-locations                  ;; Vectorof SDR indexed by LocationSDRX
      (build-vector num-locations (lambda (i) (vector-ref location-pool (+ i num-locations)))))
    (create-random-experiences         ;; Nat Nat Nat -> Experiences
      (lambda (num-exps num-sensations num-locations)
      ;; produce Object experiences, or Sequence experiences when num-locations is zero
      (let ((location-indexes (if (zero? num-locations) #f
                (build-vector (max num-locations num-sensations)
                  (if (>= num-locations num-sensations)
                    id
                    (lambda _ (random num-locations)))))))
        (build-vector num-exps (lambda _
          ;; locations are distinct but features can be repeated
          (let ((locations (if (positive? num-locations)
                               (vector-sample location-indexes num-sensations)
                               (make-vector num-sensations #f))))
            (build-vector num-sensations (lambda (sx)
                (make-sensation (random num-features) (vector-ref locations sx))))))))))
    (sequences                         ;; Vectorof Experience
      (create-random-experiences num-sequences seq-length 0))
    (objects                           ;; Vectorof Experience
      (create-random-experiences num-objects num-points num-locations))
                                                                                          ;
    ;; Stats is (ExpVectorOf (SensVectorOf (KeyListOf (CCVectorOf Number|{Nat}))))
    (stats  (build-vector (max num-objects num-sequences) (lambda _ 
                (build-vector (if (eq? figure 'f6)  (* 10 num-points)
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
              experience learn seq-ord . report)
      ;; feed each sensation of experience, report after each sensation
      (let* ( 
          (object?    (sensation-location (vector-ref experience 0)))
          (ns         (vector-length experience))
          ;; shuffle object sensations except on final training presentation
          (sensations (if (and object? (null? report))
                          (vector-sample experience ns)
                          experience))
          ;; for each sensation of object, ccs get different selection of feature+location
          ;; for sequences, all ccs get same feature, different random location
          (selects    (build-vector ns (lambda (sx)
                          (if object?
                            (build-vector num-cortical-columns (lambda (ccx)
                                (modulo (+ sx ccx) ns)))
                            (make-vector num-cortical-columns sx)))))
          (seq-locations
                  (build-vector num-cortical-columns (lambda _
                      (vector-ref (if random-seq-location  random-locations  feature-locations) 
                                  (random num-locations))))))
        (do ((sx 0 (add1 sx))) ((= sx ns))
          (let* ( 
              (select (vector-ref selects sx))
              (features
                (build-vector num-cortical-columns (lambda (ccx)
                  (if object?
                    (let ((object-sdr (vector-ref feature-pool
                            (sensation-feature (vector-ref sensations (vector-ref select ccx))))))
                      (if (and superimpose-sequence seq-ord)
                        (unique! fx=? (union1d object-sdr 
                                 (vector-ref feature-pool 
                                    (sensation-feature (vector-ref (vector-ref sequences seq-ord) sx)))))
                        object-sdr))
                    (vector-ref feature-pool 
                                (sensation-feature (vector-ref sensations sx)))))))
              (locations
                (if object?
                  (vector-refs feature-locations
                    (build-vector num-cortical-columns (lambda (ccx)
                      (sensation-location (vector-ref sensations (vector-ref select ccx))))))
                  (if location-per-sequence  seq-locations
                      (build-vector num-cortical-columns (lambda _
                          (vector-ref (if random-seq-location  random-locations  feature-locations) 
                                      (random num-locations))))))))
            (l2l4:compute patch features locations (or learn online-learning))
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
              [(l4lnp) (lengths    (lambda (l) (l2l4:get-predicted-cells        l (l2l4:l4-pop p4))) l2l4:patch-l4s)]
              [(l4lpa) (lengths    (lambda (l) (l2l4:get-predicted-active-cells l (l2l4:l4-pop p4))) l2l4:patch-l4s)]
              [(l4ac)  (vector-map (lambda (l) (l2l4:get-active-cells           l (l2l4:l4-pop p4))) (l2l4:patch-l4s patch))]
              [(l4pa)  (vector-map (lambda (l) (l2l4:get-predicted-active-cells l (l2l4:l4-pop p4))) (l2l4:patch-l4s patch))]
              [(TMlnp) (lengths    (lambda (l) (l2l4:get-predicted-cells        l (l2l4:l4-pop ss4l4))) l2l4:patch-l4s)]
              [(TMlpa) (lengths    (lambda (l) (l2l4:get-predicted-active-cells l (l2l4:l4-pop ss4l4))) l2l4:patch-l4s)]
              [(TMac)  (vector-map (lambda (l) (l2l4:get-active-cells           l (l2l4:l4-pop ss4l4))) (l2l4:patch-l4s patch))]
              [(TMpa)  (vector-map (lambda (l) (l2l4:get-predicted-active-cells l (l2l4:l4-pop ss4l4))) (l2l4:patch-l4s patch))]
              [(TXlnp) (lengths    (lambda (l) (l2l4:get-predicted-cells        l (l2l4:l4-pop ss4l23))) l2l4:patch-l4s)]
              [(TXlpa) (lengths    (lambda (l) (l2l4:get-predicted-active-cells l (l2l4:l4-pop ss4l23))) l2l4:patch-l4s)]
              [(TXac)  (vector-map (lambda (l) (l2l4:get-active-cells           l (l2l4:l4-pop ss4l23))) (l2l4:patch-l4s patch))]
              [(TXpa)  (vector-map (lambda (l) (l2l4:get-predicted-active-cells l (l2l4:l4-pop ss4l23))) (l2l4:patch-l4s patch))]
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
                  (l2l4:reset-seq patch))))
            (let ((seq (if superimpose-sequence
                          (vector-ref sequence-order (+ (* ex num-repetitions) (- num-repetitions 1)))
                          #f)))
              (have-sensations-of e #t seq ;; last training repetition: save stats
                (lambda (sx features)
                  (vector-set! (vector-ref stats ex) sx 
                    (append (get-stats-for train-keys) (list (cons 'feat features)) '())))))
            (when (positive? intersperse-noise)
              (let ((noise (vector-ref (create-random-experiences 1 seq-length 0) 0)))
                (do ((_ 0 (add1 _))) ((= _ intersperse-noise))
                  (have-sensations-of noise #t #f))))
            (l2l4:reset patch)))))     ;; reset after each object
                                                                                            ;
#;> (define (train-all-interleaved)      ;; ->
      ;; train objects & sequences, interleaving in (sqrt num-repetitions) batches
      ;; reset after inner batch of each experience
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
                    (have-sensations-of (vector-ref sequences ex) #t #f)
                    #;(l2l4:reset-seq patch))
                  (have-sensations-of (vector-ref objects ex) #t seq)))
              (when (positive? intersperse-noise)
                (let ((noise (if (zero? (random 2))
                          (create-random-experiences intersperse-noise seq-length 0)
                          (create-random-experiences intersperse-noise num-points num-locations))))
                  (do ((i 0 (add1 i))) ((= i intersperse-noise))
                    (have-sensations-of (vector-ref noise i) #t #f))))
              (l2l4:reset patch))
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
          (l2l4:reset patch))))
                                                                                            ;
#;> (define (infer-switching)            ;; ->
      ;; test mix of objects & sequences
      (let ((objectxs (vector-sample (build-vector num-objects id) 10)))
        (for-each
          (lambda (item-type ix)
            (have-sensations-of
              (vector-ref
                (if (eq? item-type 'seq) sequences  objects) (vector-ref objectxs ix))
              #f #f
              (lambda (sx features)
                (vector-set! (vector-ref stats 0) (+ (* 10 ix) sx)
                  (get-stats-for test-keys)))))
          '(seq obj seq obj seq seq obj seq obj obj) (build-list 10 id))))
                                                                                            ;
;; run experiment: 1 train network
    (if interleave-training
        (train-all-interleaved)
      (begin
        (train-all objects)
        (unless superimpose-sequence
          (train-all sequences))
        #;(train-all objects)))
;;                 2 reset and run inference
    (l2l4:reset patch)
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
            (let ((k2  (substring (symbol->string key) 0 2)))
              (if (or (string=? k2 "l4") (string=? k2 "TM") (string=? k2 "TX"))
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
                                [(l4pa)  (add-overlap-with 'l4ac)]
                                [(TMpa)  (add-overlap-with 'TMac)]
                                [(TXpa)  (add-overlap-with 'TXac)]
                                [(TMlnp TMlpa TXlnp TXlpa l4lnp l4lpa)
                                  (+ ac (vector-average (cdr stats-for-key)))]
                                [else ac])
                              ac)))
                        accs 
                        sens-for-exp)))        ;; -> (vectorof Number)
                  (build-vector (vector-length (vector-ref stats 0))
                    (lambda _ 0.0))
                  stats)) ) ] ]
                (list (symbol->string key) (list)) )))
            test-keys) ] ))
            'truncate ))
      (when display-timing
        (let ((tenths (div (+ (- (cpu-time) start-time) 50) 100))
              (sum-l2  (lambda (f)
                (vector-fold-left (lambda (sum layer)
                    (+ sum (f layer)))
                  0
                  (l2l4:patch-l2s patch))))
              (sum-l4 (lambda (f pop)
                (vector-fold-left (lambda (sum layer)
                    (+ sum (f layer pop)))
                  0
                  (l2l4:patch-l4s patch)))))
          (for-each display `(
              ,(symbol->string figure) " " ,@run-args #\newline
              ,(div tenths 10) "." ,(mod tenths 10) "s thread time\n"
              ,(sum-l4 l2l4:number-of-basal-synapses (l2l4:l4-pop ss4l4 )) "/"
              ,(sum-l4 l2l4:number-of-basal-segments (l2l4:l4-pop ss4l4 )) " ss4L4 synapses\n"
              ,(sum-l4 l2l4:number-of-basal-synapses (l2l4:l4-pop ss4l23)) "/"
              ,(sum-l4 l2l4:number-of-basal-segments (l2l4:l4-pop ss4l23)) " ss4L23 synapses\n"
              ,(sum-l4 l2l4:number-of-basal-synapses (l2l4:l4-pop p4    )) "/"
              ,(sum-l4 l2l4:number-of-basal-segments (l2l4:l4-pop p4    )) " p4 basal synapses\n"
              ,(sum-l4 l2l4:number-of-apical-synapses (l2l4:l4-pop p4   )) "/"
              ,(sum-l4 l2l4:number-of-apical-segments (l2l4:l4-pop p4   )) " p4 apical synapses\n"
              ,(sum-l2 l2l4:get-n-l2-proximal-syns) "/"
              ,(sum-l2 l2l4:get-n-l2-proximal-segs) " p23 proximal synapses\n"
              ,(sum-l2 l2l4:get-n-l2-distal-syns)   "/"
              ,(sum-l2 l2l4:get-n-l2-distal-segs)   " p23 distal synapses\n"
              ,(sum-l2 l2l4:get-n-l2-lateral-syns)  "/" 
              ,(sum-l2 l2l4:get-n-l2-lateral-segs)  " p23 lateral synapses\n"
              ))))
      ))

(define exp4a '(
    [figure         .  f4a]
    [num-sequences  .    5]
    [num-objects    .    0]
    [num-features   .   10]
    [num-locations  .  100]
    [num-repetitions .   20]
    [train-keys     . (l4ac TMac TXac)]
    [test-keys      . (l4lnp l4pa TMlnp TMpa TXlnp TXpa)]))
                                                                                            ;
(define exp5a '(
    [figure         .  f5a]
    [num-objects    .   50]
    [num-sequences  .    0]
    [num-features   .  100]
    [num-locations  .  100]
    [train-keys     . (l4ac TMac TXac)]
    [test-keys      . (l4lnp l4pa TMlnp TMpa TXlnp TXpa)]))
                                                                                            ;
(define exp6 '(
    [figure          .  f6]
    [num-sequences   .  50]
    [num-objects     .  50]
    [num-features    .  50]
    [num-locations   . 100]
    [num-repetitions .  30]
    ;[train-keys      . (l4ac TMac TXac)]
    [test-keys       . (l4lpa TMlpa TXlpa)]))
                                                                                            ;
(define run                              ;; String [ {KWarg} ] ->
  ;; can be used in repl eg (run "4a" '( [column-count . 150] ) )
  (case-lambda 
    [ (figure) (run figure '()) ]
    [ (figure options)
        (let ((start (statistics)))
          (experiment
            (case figure
              [("4a")    exp4a ]
              [("5a")    exp5a ]
              [("6")     exp6  ])
            options)
          (sstats-print (sstats-difference (statistics) start))) ] ))
  
(define (option name parameter)          ;; String [Number] -> KWarg
  ;; accept a few run options on command line
  (case name
    [("-cc")    `[column-count          . ,(string->number parameter)] ]
    [("-in")    `[intersperse-noise     . ,(string->number parameter)] ]
    [("-it")    `[interleave-training   . ,(string=? "#t"  parameter)] ]
    [("-lps")   `[location-per-sequence . ,(string=? "#t"  parameter)] ]
    [("-ncc")   `[num-cortical-columns  . ,(string->number parameter)] ]
    [("-nf")    `[num-features          . ,(string->number parameter)] ]
    [("-nl")    `[num-locations         . ,(string->number parameter)] ]
    [("-no")    `[num-objects           . ,(string->number parameter)] ]
    [("-nr")    `[num-repetitions       . ,(string->number parameter)] ]
    [("-ns")    `[num-sequences         . ,(string->number parameter)] ]
    [("-ol")    `[online-learning       . ,(string=? "#t"  parameter)] ]
    [("-rsl")   `[random-seq-location   . ,(string=? "#t"  parameter)] ]
    [("-ss")    `[superimpose-sequence  . ,(string=? "#t"  parameter)] ]))
  
(define (main command-line)              ;; {String} ->
  ;; eg $ scheme --program untangling_sequences.wp 6 -cc 150 -it
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
