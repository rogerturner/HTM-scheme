#!chezscheme

;; === HTM-scheme Untangling Sequences project  (C) 2019-2021 Roger Turner. ===
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

(library (untangling_sequences-lib)

(export experiment)

(import
          (chezscheme)
          (parameters)
          (HTM-scheme HTM-scheme algorithms htm_prelude)
          (HTM-scheme HTM-scheme algorithms htm_concept)
  (prefix (HTM-scheme HTM-scheme algorithms l2_l4_patch)  l2l4:))
                                                                                            ;
(implicit-exports #f)

;; === Types (see htm_concept.ss) ===
;; SDR            = Listof Nat
;; FeatureX       = Range num-features
;; LocationX      = Range num-locations
;; Sensation      = (FeatureX, LocationX | #f)
;; Object         = Vector CCX->Sensation
;; Sequence       = Vector CCX->Sensation
;; Experience     = Object | Sequence

(define-record-type sensation (fields
  feature                                ;; SDRX
  location))                             ;; SDRX | #f

(define (experiment exp-args run-args)   ;; KWargs KWargs ->
  ;; train l2l4 patch on object/sequence experiences, test, summarize response
  (let* (
    (args (append run-args exp-args))    ;; run-args override exp-args
    (get  (lambda (key default)
            (let ((specified (assoc key args)))
              (if specified (cdr specified) default))))
                                                                                            ;
;; model and learning parameter defaults (overridable by experiment/run parameters)
    (num-cortical-columns         (get 'num-cortical-columns  1))
    (num-minicolumns              (get 'column-count  minicolumns/macrocolumn))
    (ss4l4-cells/mcol             (get 'ss4l4-cells/mcol     15))
    (ss4l23-cells/mcol            (get 'ss4l23-cells/mcol     0))
    (p4-cells/mcol                (get 'p4-cells/mcol        15))
    (num-input-bits               (get 'num-input-bits (fxmax 5 (int<- (* num-minicolumns 0.02)))))
    (initial-permanence           (get 'initial-permanence   (perm 0.41)))
    (connected-permanence         (get 'connected-permanence (perm 0.60)))
    (permanence-increment         (get 'permanence-increment (perm 0.1 )))
    (permanence-decrement         (get 'permanence-decrement (perm 0.02)))
    (activation-threshold         (get 'activation-threshold (int<- (* num-input-bits 0.8))))
    (min-threshold                (get 'min-threshold           (int<- (* activation-threshold 0.8))))
    (reduced-basal-threshold      (get 'reduced-basal-threshold (int<- (* activation-threshold 0.8))))
    (predicted-segment-decrement  (get 'predicted-segment-decrement (perm 0.001)))
    (online-learning              (get 'online-learning      #f))
    (use-bursting-columns         (get 'use-bursting-columns #t))
                                                                                            ;
;; setup network: override algorithm default parameters
;; (scale to num-minicolumns - works for 100-400 minicolumns?)
    (l2-cells/mcol                10)
    (l2-sdr-size                  (int<- (* 1.5 num-input-bits)))
    (l6-cells/mcol                5)
    (external-input-size          (fxmin 2000 (* num-minicolumns l6-cells/mcol)))  ;; p6l4
    (location-bits                (int<- (* external-input-size 0.01)))
    (l2-sample-size-proximal      (int<- (* 1.0 l2-sdr-size)))
    (l2-sample-size-distal        (* 7 l2-sample-size-proximal))
    (l2-overrides  (get 'l2-overrides `(
      [initial-proximal-permanence    . ,initial-permanence]
      [connected-permanence-proximal  . ,connected-permanence]
      [initial-distal-permanence      . ,initial-permanence]
      [connected-permanence-distal    . ,connected-permanence]
      [sample-size-proximal           . ,l2-sample-size-proximal]
      [min-threshold-proximal         . ,(int<- (* 0.8 num-input-bits))]
      [predicted-inhibition-threshold . ,l2-sample-size-proximal]
      [sample-size-distal             . ,l2-sample-size-distal]
      [activation-threshold-distal    . ,(int<- (* 0.8 num-input-bits))]
      [sdr-size                       . ,l2-sdr-size]
      [min-sdr-size                   . ,(int<- (* 0.75 l2-sdr-size))]
      [online-learning                . ,online-learning])))
    (l4-overrides  (get 'l4-overrides `(
      [perm-trim-threshold                . ,(perm 0.31)]
      [use-bursting-columns-input         . ,use-bursting-columns]
      [sample-size                        . ,(int<- (* 1.5 num-input-bits))]
      [initial-permanence                 . ,initial-permanence]
      [connected-permanence               . ,connected-permanence]
      [permanence-increment               . ,permanence-increment]
      [permanence-decrement               . ,permanence-decrement]
      [activation-threshold               . ,activation-threshold]
      [min-threshold                      . ,min-threshold]
      [reduced-basal-threshold            . ,reduced-basal-threshold]
      [apical-predicted-segment-decrement . ,predicted-segment-decrement]
      [basal-predicted-segment-decrement  . ,predicted-segment-decrement])))
    (patch 
      (l2l4:make-patch num-cortical-columns num-minicolumns l2-cells/mcol
        ss4l4-cells/mcol ss4l23-cells/mcol p4-cells/mcol l6-cells/mcol
        l2-overrides l4-overrides))
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
    (random-sdr (lambda (size bits)      ;; Nat Nat -> ( -> SDR)
        ;; produce function to generate SDR of size with sparsity bits/size
        (lambda ()
          (sort fx<? (u32-sample (iota size) bits)))))
    (build-distinct (lambda (n b ->sdr)  ;; Nat Nat (-> SDR) -> (Vector N->SDR)
        ;; produce n vector of b on-bit SDRs with minimum overlaps (no more than b/5 bits)
        (let ([v (make-vector n)])
          (do ([vx 0 (fx1+ vx)]) ((fx=? vx n) v)
            (let try1 ([overlap 0])
              (if (fx>? overlap 1 #;(fx1+ (fxdiv b 5)))  (error #f "can't generate distinct features or locations")
                (let try2 ([tries 0])
                  (if (fx>? tries 1000) (try1 (fx1+ overlap))
                    (let ([sdr (->sdr)])
                      (do ([vx- 0 (fx1+ vx-)])
                          ((or (fx=? vx- vx)
                               (fx>? (length (intersect1d (vector-ref v vx-) sdr)) overlap))
                            (if (fx=? vx- vx)  (vector-set! v vx sdr)
                              (try2 (fx1+ tries))))))))))))))
    (feature-pool                        ;; Vector FeatureX->SDR
      (build-distinct num-features num-input-bits (random-sdr num-minicolumns num-input-bits)))
    (location-pool                       ;; Vector LocationX->SDR
      (build-distinct (* 2 num-locations) location-bits (random-sdr external-input-size location-bits)))
    (random-locations                       ;; Vector LocationX->SDR
      (build-distinct num-locations location-bits (random-sdr external-input-size location-bits)))
    (create-random-experiences         ;; Nat Nat Nat -> Experiences
      (lambda (num-exps num-sensations num-locations)
        ;; produce Object experiences, or Sequence experiences when num-locations is zero
        (build-vector num-exps (lambda _
          ;; locations are distinct but features can be repeated
          (let ([locations (if (positive? num-locations)
                             (vector-sample (indexes feature-pool) num-sensations)
                             (make-vector num-sensations #f))]
                [features  (build-vector num-sensations (lambda _ (random num-features)))])
            (build-vector num-sensations (lambda (sx)
                (make-sensation (vector-ref features sx) (vector-ref locations sx)))))))))
    (sequences                         ;; Vector SequenceX->Experience
      (create-random-experiences num-sequences seq-length 0))
    (objects                           ;; Vector ObjectX->Experience
      (create-random-experiences num-objects num-points num-locations))
                                                                                          ;
    ;; Stats is (Vector ExperienceX->(Vector SensationX->(KeyListOf (Vector CCX->Number|{Nat}))))
    (stats  (build-vector (max num-objects num-sequences) (lambda _ 
                (build-vector (if (eq? figure 'f6)  (* 10 num-points)
                                  (max num-points seq-length))
                              (lambda _ (list))))))
    (display-timing (get 'display-timing  #t))
    (train-keys     (get 'train-keys     '()))
    (test-keys      (get 'test-keys      '()))
    (test-keys      (filter (lambda (key)
                        (let ([k2 (substring (symbol->string key) 0 2)])
                          (case k2
                            [("TM") (positive? ss4l4-cells/mcol)]
                            [("TX") (positive? ss4l23-cells/mcol)]
                            [("l4") (positive? p4-cells/mcol)]
                            [else #t])))
                      test-keys))
    (f6-objectxs    #f)
    (start-time     (cpu-time))
    )

#;> (define (have-sensations-of          ;; Experience Bool ( -> (SDR Nat -> SDR)) [(Nat -> )] ->
              experience learn . report)
      ;; feed each sensation of experience to model, optional report after each sensation
      (let* ( 
          [object?    (sensation-location (vector-ref experience 0))]
          [ns         (vector-length experience)]
          ;; shuffle object sensations except on final training presentation
          [sensations (if (and object? (null? report))
                          (vector-sample experience ns)
                          experience)]
          ;; for each sensation of object, ccs get different selection of feature+location
          ;; for sequences, all ccs get same feature, different random location
          [sx->>sx    (build-vector ns (lambda (sx)    ;; Vector SensationX->(Vector CCX->SensationX)
                          (if object?
                            (build-vector num-cortical-columns (lambda (ccx)
                                (modulo (fx+ (fxdiv ccx (ceiling (/ num-cortical-columns ns))) sx) ns)))
                            (make-vector num-cortical-columns sx))))]
          [seq-locations
                  (build-vector num-cortical-columns (lambda _
                      (vector-ref location-pool #;(if random-seq-location  random-locations  object-locations) 
                                  (+ num-locations (random num-locations)))))])
        (do ([sx 0 (fx1+ sx)]) ((fx=? sx ns))
          (let* ( 
              [ccx->sx (vector-ref sx->>sx sx)]
              [features
                (build-vector num-cortical-columns (lambda (ccx)
                  (if object?
                    (vector-ref feature-pool
                      (sensation-feature (vector-ref sensations (vector-ref ccx->sx ccx))))
                    (vector-ref feature-pool 
                      (sensation-feature (vector-ref sensations sx))))))]
              [locations
                #;(if object?
                  (vector-refs location-pool #;object-locations
                    (build-vector num-cortical-columns (lambda (ccx)
                      (sensation-location (vector-ref sensations (vector-ref ccx->sx ccx))))))
                  (if location-per-sequence  seq-locations
                      (build-vector num-cortical-columns (lambda _
                          (vector-ref location-pool #;(if random-seq-location  random-locations  object-locations) 
                                    (+ num-locations (random num-locations)))))))
                (build-vector num-cortical-columns (lambda (ccx)
                  (if object?
                    (vector-ref location-pool #;object-locations
                      (sensation-location (vector-ref sensations (vector-ref ccx->sx ccx))))
                    (vector-ref location-pool #;(if random-seq-location  random-locations  object-locations) 
                        (+ num-locations (random num-locations))))))])
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
              [(l4lnp) (lengths    (lambda (l) (l2l4:get-predicted-cells        l 'p4    )) l2l4:patch-l4s)]
              [(l4lpa) (lengths    (lambda (l) (l2l4:get-predicted-active-cells l 'p4    )) l2l4:patch-l4s)]
              [(l4ac)  (vector-map (lambda (l) (l2l4:get-active-cells           l 'p4    )) (l2l4:patch-l4s patch))]
              [(l4pa)  (vector-map (lambda (l) (l2l4:get-predicted-active-cells l 'p4    )) (l2l4:patch-l4s patch))]
 ;              [(l3lpa) (lengths    (lambda (l) (l2l4:get-predicted-active-cells-l3 l)) l2l4:patch-l3s)]
 ;              [(l3ac)  (vector-map (lambda (l) (l2l4:get-active-cells-l3        l))                  (l2l4:patch-l3s patch))]
 ;              [(l3pa)  (vector-map (lambda (l) (l2l4:get-predicted-active-cells-l3 l))               (l2l4:patch-l3s patch))]
              [(TMlnp) (lengths    (lambda (l) (l2l4:get-predicted-cells        l 'ss4l4 )) l2l4:patch-l4s)]
              [(TMlpa) (lengths    (lambda (l) (l2l4:get-predicted-active-cells l 'ss4l4 )) l2l4:patch-l4s)]
              [(TMac)  (vector-map (lambda (l) (l2l4:get-active-cells           l 'ss4l4 )) (l2l4:patch-l4s patch))]
              [(TMpa)  (vector-map (lambda (l) (l2l4:get-predicted-active-cells l 'ss4l4 )) (l2l4:patch-l4s patch))]
              [(TXlnp) (lengths    (lambda (l) (l2l4:get-predicted-cells        l 'ss4l23)) l2l4:patch-l4s)]
              [(TXlpa) (lengths    (lambda (l) (l2l4:get-predicted-active-cells l 'ss4l23)) l2l4:patch-l4s)]
              [(TXac)  (vector-map (lambda (l) (l2l4:get-active-cells           l 'ss4l23)) (l2l4:patch-l4s patch))]
              [(TXpa)  (vector-map (lambda (l) (l2l4:get-predicted-active-cells l 'ss4l23)) (l2l4:patch-l4s patch))]
              )))
        keys))
                                                                                            ;
#;> (define (train-all es)               ;; Experiences ->
      ;; train on each e num-repetitions times
      (do-with-progress (vector-length es)
        (lambda (ex)
          (let ((e (vector-ref es ex)))
            (do ([r 0 (fx1+ r)]) ((fx>=? r (fx1- num-repetitions)))
              (have-sensations-of e #t)
              (unless (sensation-location (vector-ref e 0))
                (l2l4:reset-seq patch))) ;; reset TM between presentations of sequence
            (when (positive? intersperse-noise)
              (let ((noise (vector-ref (create-random-experiences 1 seq-length 0) 0)))
                (do ((_ 0 (fx1+ _))) ((fx=? _ intersperse-noise))
                  (have-sensations-of noise #t))))
            (l2l4:reset patch)))))       ;; reset after each object
                                                                                            ;
#;> #;(define (train-all-interleaved)      ;; ->
      ;; train objects & sequences
      ;; reset after inner batch of each experience
      (assert (= num-objects num-sequences))
      (let ([seq-tries (make-vector num-sequences num-repetitions)]
            [obj-tries (make-vector num-objects   num-repetitions)])
        (let next-pass ([passes-left 8])
          (when (and (positive? passes-left)
                     (exists (lambda (st ot)
                         (or (> st 4) (> ot 4)))
                       (vector->list seq-tries) (vector->list obj-tries)))
            (do ([ex 0 (+ ex 1)]) ((= ex num-objects))
              ;(newline)
              (let (
                [train-one
                  (lambda (seq/obj tries pop)
                    (when (= (vector-ref tries ex) num-repetitions)
                      (l2l4:reset patch))
                    (have-sensations-of (vector-ref seq/obj ex) #t #f)
                    (have-sensations-of (vector-ref seq/obj ex) #t #f)
                    (have-sensations-of (vector-ref seq/obj ex) #t #f)
                    (let train-one ([tx (quotient (vector-ref tries ex) 4)])
                      (cond
                        [ (zero? tx) ]
                        [ (for-all (lambda (l)  ;; all ccs
                              (= num-input-bits (length (l2l4:get-active-cells l pop))))
                            (vector->list (l2l4:patch-l4s patch)))
                          (vector-set! tries ex 0) ]
                        [ else 
                          (have-sensations-of (vector-ref seq/obj ex) #t #f)
                          (vector-set! tries ex (- (vector-ref tries ex) 1))
                          (train-one (- tx 1)) ])))])
                (train-one sequences seq-tries 'ss4l4)
                (train-one objects   obj-tries 'p4) ))
            (display seq-tries) (display obj-tries)
            (next-pass (- passes-left 1))))))
                                                                                            ;

(define nt 0)

#;> #;(define (train-all-interleaved)      ;; ->
      ;; train objects & sequences
      ;; reset only before first encounter of each experience
      (assert (= num-objects num-sequences))
      (let ([seq-tries (make-vector num-sequences num-repetitions)]
            [obj-tries (make-vector num-objects   num-repetitions)])
        (let next-pass ([passes-left num-repetitions])
          (when (and (positive? passes-left)
                     (exists (lambda (st ot)
                         (or (> st 0) (> ot 0)))
                       (vector->list seq-tries) (vector->list obj-tries)))
            (do ([ex 0 (+ ex 1)]) ((= ex num-objects))
              (vector-set! seq-tries ex (max passes-left (vector-ref seq-tries ex)))
              (vector-set! obj-tries ex (max passes-left (vector-ref obj-tries ex))))
              
            (do ([ex 0 (+ ex 1)]) ((= ex num-objects))
              (let (
                [train-one-item
                  (lambda (seq/obj tries pop)
                    (when (= (vector-ref tries ex) num-repetitions)
                      (l2l4:reset-l2 patch))
                    (have-sensations-of (vector-ref seq/obj ex) #t #f)
                    (vector-set! tries ex (- (vector-ref tries ex) 1))
                    (set! nt (+ nt 1))
                    (let train-until ([tx (vector-ref tries ex)])
                      (cond
                        [ (<= tx 0) ]
                        [ (for-all (lambda (l)  ;; all ccs
                              (= num-input-bits (length (l2l4:get-active-cells l pop))))
                            (vector->list (l2l4:patch-l4s patch))) ]
                        [ else 
                          (have-sensations-of (vector-ref seq/obj ex) #t #f)
                          (vector-set! tries ex (- (vector-ref tries ex) 1))
                          (set! nt (+ nt 1))
                          (train-until (- tx 1)) ])))])
                (train-one-item sequences seq-tries 'ss4l4)
                (train-one-item objects   obj-tries (if (positive? p4-cells/mcol) 'p4 'ss4l23))
                 ))
            (next-pass (- passes-left 1))))
            (display seq-tries) (display obj-tries) (newline)
            (display nt) (display " trains ")
            
            ))
                                                                                            ;
#;> (define (train-all-interleaved)      ;; ->
      ;; train objects & sequences
      ;; reset l2 only before first encounter of each experience
      ;; reset sequence memory (ss4l4) before each sequence
      (assert (= num-objects num-sequences))
    (let ()
      (define (train-all-items max-tries)
        (define (train-one-item seq/obj ex tries pop)
          ;; train-one-item
          (when (= tries num-repetitions)
            (l2l4:reset-l2 patch))
          (let train-up-to ([tx tries])
            (when (> tx 8)
              (have-sensations-of (vector-ref seq/obj ex) #t)
              (set! nt (+ nt 1)))
            (unless (or (zero? tx)
                (for-all (lambda (l)  ;; all ccs
                    (= num-input-bits (length (l2l4:get-active-cells l pop))))
                  (vector->list (l2l4:patch-l4s patch))))
              (have-sensations-of (vector-ref seq/obj ex) #t)
              (set! nt (+ nt 1))
              (train-up-to (- tx 1)))))
        ;; train-all-items
        (do ([ex 0 (+ ex 1)]) ((= ex num-objects))
          (l2l4:reset-l4 patch)
          (train-one-item sequences ex max-tries 'ss4l4)
          (l2l4:reset-l4 patch)
          (train-one-item objects   ex max-tries (if (positive? p4-cells/mcol) 'p4 'ss4l23))))
      (do ([max-tries num-repetitions
                      (if (> max-tries 16) (quotient max-tries 2) (- max-tries 1))])
          ((zero? max-tries))
        (train-all-items max-tries))
      (train-all-items 1)
      (display nt) (display " trains ") (newline)))
                                                                                            ;
#;> (define (infer-all es)               ;; Experiences ->
      ;; test each experience and save stats for each sensation
      (do-with-progress (vector-length es)
        (lambda (ex)
          (have-sensations-of (vector-ref es ex) #f
            (lambda (sx features)
              (vector-set! (vector-ref stats ex) sx
                (append (get-stats-for test-keys)
                        (vector-ref (vector-ref stats ex) sx)
                        '()))))
          (l2l4:reset patch))))
                                                                                            ;
#;> (define (infer-switching)            ;; ->
      ;; test mix of objects & sequences
      (for-each
        (lambda (item-type ix)
          (have-sensations-of
            (vector-ref
              (if (eq? item-type 'seq) sequences  objects) (vector-ref f6-objectxs ix))
            #f
            (lambda (sx features)
              (vector-set! (vector-ref stats 0) (fx+ (fx* num-points ix) sx)
                (append (get-stats-for test-keys)
                        (vector-ref (vector-ref stats 0) (fx+ (fx* num-points ix) sx))
                        '()))))
          #;(l2l4:reset patch))
        '(seq obj seq obj seq seq obj seq obj obj) (iota 10)))
                                                                                            ;
#;> (define (stats->plot ccx)
      [cons [list "figure" (symbol->string figure)]
      [cons [list "using"  (cons `[num-minicolumns . ,minicolumns/macrocolumn] run-args)]
      (map (lambda (key)             ;; each stats key (Symbol -> )
          (let ([k2  (substring (symbol->string key) 0 2)])
            (if (or (string=? k2 "l3") (string=? k2 "l4") (string=? k2 "TM") (string=? k2 "TX"))
              [cons (symbol->string key)
              ;; stats is (ExpVectorOf (SensVectorOf (KeyListOf (CCVectorOf Number))))
              [list (vector->list
                (vector-map (lambda (x)
                    (if (eq? figure 'f6)  x   ;; average of exps (for 4a/5a)
                      (int<- (/ x (fx+ num-objects num-sequences)))))
                  (vector-fold-left (lambda (accs sens-for-exp)  ;; (listof Number) (SensListOf (KeyListOf (vectorof Number)))
                      ;; ..accumulated objects|sequences
                      (if (null? (vector-ref sens-for-exp 0))  accs
                        (vector-map-x (lambda (ac keys-stats sensx)    ;; ..each sensation, average ccs
                            ;; Number (Alistof (vectorof Number)) -> Number
                            (let ([stats-for-key (assq key keys-stats)])
                              (define (overlap t a)
                                (+ (length (intersect1d t a)) (length (setdiff1d t a)) #;(length (setdiff1d a t)))
                                #;(+ (length (intersect1d t a)) (* -10 (length (setdiff1d t a))) (* -10 (length (setdiff1d a t)))))
                              (define (overlap-with train-key)
                                (let ([trained (assq train-key keys-stats)])
                                  (if trained
                                    #;(vector-map (lambda (t a)
                                        (length a)
                                        #;(+ (length (intersect1d t a))
                                            (if #f #;(zero? (modulo sensx 10))  0
                                                (+ #;(length (setdiff1d t a))
                                                   (length (setdiff1d a t))))))
                                      (cdr trained) (cdr stats-for-key))
                                    #;(vector (length (vector-ref (cdr stats-for-key) ccx)))
                                    (vector (if (and (zero? (modulo sensx 10))
                                                     (memv (quotient sensx 10) '(0 2 4 5 7)))
                                              0
                                              (overlap (vector-ref (cdr trained) ccx)
                                                       (vector-ref (cdr stats-for-key) ccx))))
                                    ac)))
                              (define (vector-binarize v)
                                (vector-map (lambda (e)
                                    (cond
                                      [(fx<=? (fx1- num-input-bits) e) (fxmax e num-input-bits)]
                                      [(fx<=? e 1) 0]
                                      [else e]))
                                  v))
                              (if stats-for-key
                                (case key
                                  [(l3pa)  (overlap-with 'l3ac)]
                                  [(l4pa)  (overlap-with 'l4ac)]
                                  [(TMpa)  (overlap-with 'TMac)]
                                  [(TXpa)  (overlap-with 'TXac)]
                                  [(TMlnp TMlpa TXlnp TXlpa l4lnp l4lpa l3lpa)
                                    (cdr stats-for-key)
                                    #;(vector-binarize (cdr stats-for-key))
                                    #;(+ ac (vector-average (cdr stats-for-key)))]
                                  [else ac])
                                ac)))
                          accs 
                          sens-for-exp)))        ;; -> (vectorof Number)
                    (build-vector (vector-length (vector-ref stats 0))
                        (lambda _ 0.0))
                    stats)) ) ] ]
              (list (symbol->string key) (list)) )))
        test-keys) ]] )
                                                                                            ;
;; run experiment: 1 train network
    (if interleave-training
        (train-all-interleaved)
      (begin
        (train-all objects)
        (unless superimpose-sequence
          (train-all sequences))))
;;                 2 save objects/sequences and trained stats
    (when (and (eq? figure 'f6) (positive? num-repetitions))
      (set! f6-objectxs (vector-sample (build-vector num-objects id) 10))
      (for-each
        (lambda (item-type ix)
          (have-sensations-of
            (vector-ref
              (if (eq? item-type 'seq) sequences  objects) (vector-ref f6-objectxs ix))
            #t
            (lambda (sx features)
              (vector-set! (vector-ref stats 0) (fx+ (fx* num-points ix) sx)
                (get-stats-for train-keys)))))
        '(seq obj seq obj seq seq obj seq obj obj) (iota 10)))
;;                 3 reset and run inference
    (l2l4:reset patch)
    (when (positive? num-repetitions)
      (case figure
        [(f6)
          (infer-switching) ]
        [else
          (infer-all objects)
          (infer-all sequences) ] ))
;;                 4 save statistics for plots
    (with-output-to-file "experiment.data" (lambda ()
        (write (map stats->plot (iota num-cortical-columns))))
      'truncate )
    (when display-timing
      (let ([tenths (div (+ (- (cpu-time) start-time) 50) 100)])
        (for-each display `(,(div tenths 10) "." ,(mod tenths 10) "s thread time\n")))
      (let ([t1 (cost-center-time cost-center-1)]
            [t2 (cost-center-time cost-center-2)]
            [tz (make-time 'time-duration 0 0)])
        (when (time>? t1 tz)
          (for-each display `(,(time-second t1) "." ,(div (time-nanosecond t1) 10000000) "s in cost center 1\n")))
        (when (time>? t2 tz)
          (for-each display `(,(time-second t2) "." ,(div (time-nanosecond t2) 10000000) "s in cost center 2\n"))))
      (for-each display `(,(symbol->string figure) " " ,@run-args "        \n"))
      (l2l4:print-statistics patch))
      )
  )

)