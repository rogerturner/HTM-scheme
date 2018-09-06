#!r6rs

;; === HTM-scheme Apical Tiebreak Sequence Memory Copyright 2017 Roger Turner. ===
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

  ;; Extends the HTM-scheme Apical Tiebreak Temporal Memory library.
  ;; Translated from numenta htmresearch/.../apical_tiebreak_temporal_memory.py --
  ;; see comments there for descriptions of functions and parameters.
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.

(library (HTM-scheme HTM-scheme algorithms apical_tiebreak_sequence_memory)
                                                                                            ;
(export
  tm:permanence
  tm:synapse
  tm:segment
  (rename
    (tm:get-active-cols              atsm:get-active-cols)
    (tm:get-active-cells             atsm:get-active-cells)
    (tm                              atsm:tm)
    (make-tm                         atsm:construct)
    (reset                           atsm:reset)
    (compute                         atsm:compute)
    (get-predicted-cells             atsm:get-predicted-cells)
    (get-next-predicted-cells        atsm:get-next-predicted-cells)
    (get-next-basal-predicted-cells  atsm:get-next-basal-predicted-cells)
    (get-next-apical-predicted-cells atsm:get-next-apical-predicted-cells)))
                                                                                            ;
(import 
  (rnrs)
  (HTM-scheme HTM-scheme algorithms htm_prelude)
  (HTM-scheme HTM-scheme algorithms apical_tiebreak_temporal_memory))
                                                                                            ;
(define-record-type tm                   ;; TM
  (parent tm:tm)
  (fields
    (mutable prev-apical-input)             ;; (listof CellX)
    (mutable prev-apical-growth-candidates) ;; (listof CellX)
    (mutable prev-predicted-cells))         ;; (listof CellX)
  (protocol
    (lambda (pargs->new)                 ;; Nat Nat (listof KWarg) -> TM
      (lambda (column-count cells-per-column . kwargs)
          (apply (apply pargs->new (append (list column-count cells-per-column) kwargs))
                 (key-word-args
                    kwargs
                    atsm-defaults #f))))))
                                                                                           ;
(define atsm-defaults                    ;; (listof KWarg)
  `(
    [prev-apical-input                  . ()]
    [prev-apical-growth-candidates      . ()]
    [prev-predicted-cells               . ()]))
    
                                                                                            ;
(define (reset tm)                       ;; TM ->
  (tm:reset tm)
  (tm-prev-apical-input-set! tm             '())
  (tm-prev-apical-growth-candidates-set! tm '())
  (tm-prev-predicted-cells-set! tm          '()))
                                                                                            ;
(define compute                          ;; TM {ColX} {CellX} [{CellX}] Boolean ->
  (case-lambda
  [(tm active-columns apical-input apical-growth-candidates learn)
  ;; Perform one timestep. Activate the specified columns, using the predictions
  ;; from the previous timestep, then learn. Then form a new set of predictions
  ;; using the new active cells and the apicalInput.
    (let ((apical-growth-candidates
            (if (null? apical-growth-candidates)
                apical-input
                apical-growth-candidates)))
      (tm-prev-predicted-cells-set! tm (tm-predicted-cells tm))
      (tm:activate-cells tm active-columns (tm:get-active-cells tm) (tm-prev-apical-input tm)
                      (tm:get-winner-cells tm) (tm-prev-apical-growth-candidates tm) learn)
      (tm:depolarize-cells tm (tm:get-active-cells tm) apical-input learn)
      (tm-prev-apical-input-set! tm apical-input)
      (tm-prev-apical-growth-candidates-set! tm apical-growth-candidates))]
  [(tm active-columns apical-input learn)
    (compute tm active-columns apical-input apical-input learn)]
  [(tm active-columns learn)
    (compute tm active-columns '() learn)]))
                                                                                           ;
(define (get-predicted-cells tm)         ;; TM -> {CellX}
  ;; The prediction from the previous timestep
  (tm-prev-predicted-cells tm))
                                                                                            ;
(define (get-next-predicted-cells tm)    ;; TM -> {CellX}
  ;; The prediction for the next timestep
  (tm-predicted-cells tm))
                                                                                            ;
(define (get-next-basal-predicted-cells tm)  ;; TM -> {CellX}
  ;; Cells with active basal segments
  (unique (map-segments-to-cells (tm-active-basal-segments tm))))
                                                                                            ;
(define (get-next-apical-predicted-cells tm) ;; TM -> {CellX}
  ;; Cells with active apical segments
  (unique (map-segments-to-cells (tm-active-apical-segments tm))))

;; === Smoke tests ===
                                                                                          ;
(define tests "")
(define failures #f)
                                                                                            ;
(define (test name)
  (set! tests (string-append tests "\n" name)))
                                                                                            ;
(define-syntax expect                    ;; ((X ... -> Y) X ...) Y -> [error]
  ;; check that function application(s) to arguments match expected values
  (lambda (x)                            
    (syntax-case x ()                                  ;; [expect ([fn args] expected ) ... ]
      [ (_ (expr expected) ...)                        ;; expr matches [fn args]
        #'(begin (let ((result expr))                  ;; eval expr just once
                   (if (equal? result expected)
                      (set! tests (string-append tests " ok"))
                      (begin 
                        (for-each display `(,tests #\newline "**" expr #\newline
                          "  expected: " ,expected #\newline 
                          "  returned: " ,result  #\newline))
                        (set! failures #t)
                        (set! tests "")))) ...)])))
                                                                                            ;
;(define print-fields
;  ;; Excerpt From: R. Kent Dybvig. The Scheme Programming Language, 4th Edition.
;  (lambda (r)
;    (unless (record? r)
;      (assertion-violation 'print-fields "not a record" r))
;    (let loop ([rtd (record-rtd r)])
;      (let ([prtd (record-type-parent rtd)])
;        (when prtd (loop prtd)))
;      (let* ([v (record-type-field-names rtd)]
;             [n (vector-length v)])
;        (do ([i 0 (+ i 1)])
;            ((= i n))
;          (write (vector-ref v i))
;          (display "=")
;          (write ((record-accessor rtd i) r))
;          (newline))))))
                                                                                            ;
(define (make-test-tm . options)
  (random-seed! 42)
  (apply make-tm 
    (append 
      (list 32 4)
      options
      `([activation-threshold . 3]
        [min-threshold        . 2]
        [sample-size          . 3]
        [basal-predicted-segment-decrement  . 0]
        [apical-predicted-segment-decrement . 0]
        [initial-permanence . ,(tm:permanence 0.5)] ;; so that grow-synapses creates connected synapses
      ))))
                                                                                            ;
(define (create-basal-segment tm segx)
  (create-segment tm segx (tm-basal-connections tm)))
                                                                                            ;
(define (create-basal-synapses tm seg pre-cells)
  (grow-synapses tm seg (length pre-cells) pre-cells (tm-basal-pre-index tm)))

(test "ActivateCorrectlyPredictiveCells")
  (let* ( (tm (make-test-tm))
        (active-segment (create-basal-segment tm 4)))
  (create-basal-synapses tm active-segment '(0 1 2 3))
  (compute tm '(0) #t)
  (compute tm '(1) #t)
  [expect ([get-predicted-cells tm] '(4))
          ([tm:get-active-cells    tm] '(4))])
                                                                                            ;
(test "BurstUnpredictedColumns")
  (let* ( (tm (make-test-tm)))
  (compute tm '(0) #t)
  [expect ([get-predicted-cells tm] '())
          ([tm:get-active-cells    tm] '(0 1 2 3))])
                                                                                            ;
(test "ZeroActiveColumns")
  (let* ( (tm (make-test-tm `[predicted-segment-decrement . ,(tm:permanence 0.02)]))
        (segment (create-basal-segment tm 4)))
  (create-basal-synapses tm segment '(0 1 2 3))
  (compute tm '(0) #t)
  [expect ((null? [tm:get-active-cells tm]) #f)
          ((null? [tm:get-winner-cells tm])  #f)]
  (compute tm '() #t)
  [expect ((null? [tm:get-active-cells tm]) #t)
          ((null? [tm:get-winner-cells tm])  #t)])
                                                                                            ;
(test "PredictedActiveCellsAreAlwaysWinners")
  (let* ( (tm (make-test-tm))
        (active-segment-1 (create-basal-segment tm 4))
        (active-segment-2 (create-basal-segment tm 6)))
  (create-basal-synapses tm active-segment-1 '(0 1 2))
  (create-basal-synapses tm active-segment-2 '(0 1 2))
  (compute tm '(0) '() '() #f)
  (compute tm '(1) '() '() #f)
  [expect ([tm:get-winner-cells tm]  '(4 6))])
                                                                                            ;
(test "ChooseOneWinnerCellInBurstingColumn")
  (let* ( (tm (make-test-tm)))
  (compute tm '(0) '() '() #f)
  [expect ((length [tm:get-winner-cells tm]) 1)]
  (unless (null? [tm:get-winner-cells tm])
    [expect ((not (member (car [tm:get-winner-cells tm]) '(0 1 2 3))) #f)]))
                                                                                            ;
(test "ReinforceCorrectlyActiveSegments")
  (let* ( (tm (make-test-tm
              `[sample-size                 . 4]
              `[permanence-decrement        . ,(tm:permanence 0.08)]
              `[predicted-segment-decrement . ,(tm:permanence 0.02)]))
        (active-segment (create-basal-segment tm 5)))
  (create-basal-synapses tm active-segment '(0 1 2 81))
  (compute tm '(0) #t)
  (compute tm '(1) #t)
  [expect ((vector-map syn-perm [seg-synapses active-segment]) '#(6000 6000 6000 4200))])
                                                                                            ;
(test "ReinforceSelectedMatchingSegmentInBurstingColumn")
  (let* ( (tm (make-test-tm
              `[permanence-decrement . ,(tm:permanence 0.08)]
              `[initial-permanence   . ,(tm:permanence 0.3)]))
        (selected-matching-segment (create-basal-segment tm 4))
        (other-matching-segment    (create-basal-segment tm 5)))
  (create-basal-synapses tm selected-matching-segment '(0 1 2 81))
  (create-basal-synapses tm other-matching-segment    '(0 1 81))
  (compute tm '(0) #t)
  (compute tm '(1) #t)
  [expect ((vector-map syn-perm [seg-synapses selected-matching-segment]) '#(4000 4000 4000 2200))])
                                                                                          ;
(test "NoChangeToNonselectedMatchingSegmentsInBurstingColumn")
  (let* ( (tm (make-test-tm
              `[permanence-decrement . ,(tm:permanence 0.08)]
              `[initial-permanence   . ,(tm:permanence 0.3)]))
        (selected-matching-segment (create-basal-segment tm 4))
        (other-matching-segment    (create-basal-segment tm 5)))
  (create-basal-synapses tm selected-matching-segment '(0 1 2 81))
  (create-basal-synapses tm other-matching-segment    '(0 1 81))
  (compute tm '(0) #t)
  (compute tm '(1) #t)
  [expect ((vector-map syn-perm [seg-synapses other-matching-segment]) '#(3000 3000 3000))])
                                                                                          ;
(test "NoChangeToMatchingSegmentsInPredictedActiveColumn")
  (let* ( (tm (make-test-tm))
        (active-segment (create-basal-segment tm 4))
        (matching-segment-on-same-cell  (create-basal-segment tm 4))
        (matching-segment-on-other-cell (create-basal-segment tm 5)))
  (create-basal-synapses tm active-segment '(0 1 2 3))
  (create-basal-synapses tm matching-segment-on-same-cell '(0 1))
  (seg-synapses-set! matching-segment-on-same-cell (vector (tm:synapse 0 3000) (tm:synapse 1 3000)))
  (create-basal-synapses tm matching-segment-on-other-cell '(0 1))
  (seg-synapses-set! matching-segment-on-other-cell (vector (tm:synapse 0 3000) (tm:synapse 1 3000)))
  (compute tm '(0) #t)
  (compute tm '(1) #t)
  [expect ([get-predicted-cells tm] '(4))
          ((vector-map syn-perm [seg-synapses matching-segment-on-same-cell])  '#(3000 3000))
          ((vector-map syn-perm [seg-synapses matching-segment-on-other-cell]) '#(3000 3000))])
                                                                                          ;
(test "NoNewSegmentIfNotEnoughWinnerCells")
  (let* ( (tm (make-test-tm
              `[sample-size                 . 2])))
  (compute tm '() #t)
  (compute tm '(0) #t)
  [expect ([tm-next-flatx tm] 0)])
                                                                                          ;
(test "NewSegmentAddSynapsesToSubsetOfWinnerCells")
  (let* ( (tm (make-test-tm
              `[sample-size        . 2]
              `[initial-permanence . ,(tm:permanence 0.21)])))
  (compute tm '(0 1 2) #t)
  (let ((prev-winner-cells [tm:get-winner-cells tm]))
    [expect ((length prev-winner-cells) 3)]
    (compute tm '(4) #t)
    [expect ((length [tm:get-winner-cells tm]) 1)]
    (let* ( (segments (vector-ref [tm-basal-connections tm] (car [tm:get-winner-cells tm])))
            (synapses (seg-synapses (car segments))))
      [expect ((length segments) 1)
              ((vector-length synapses) 2)]
      (vector-for-each
        (lambda (synapse)
          [expect ((syn-perm synapse) 2100)
                  ((not (member (syn-prex synapse) prev-winner-cells)) #f)])
        synapses))))
                                                                                          ;
(test "NewSegmentAddSynapsesToAllWinnerCells")
  (let* ( (tm (make-test-tm
              `[sample-size        . 4]
              `[initial-permanence . ,(tm:permanence 0.21)])))
  (compute tm '(0 1 2) #t)
  (let ((prev-winner-cells [tm:get-winner-cells tm]))
    [expect ((length prev-winner-cells) 3)]
    (compute tm '(4) #t)
    [expect ((length [tm:get-winner-cells tm]) 1)]
    (let* ( (segments (vector-ref [tm-basal-connections tm] (car [tm:get-winner-cells tm])))
            (synapses (seg-synapses (car segments))))
      [expect ((length segments) 1)
              ((vector-length synapses) 3)]
      (vector-for-each
        (lambda (synapse)
          [expect ((syn-perm synapse) 2100)
                  ((not (member (syn-prex synapse) prev-winner-cells)) #f)])
        synapses)
      [expect ((vector->list (vector-map syn-prex synapses)) prev-winner-cells)])))
                                                                                          ;
(test "MatchingSegmentAddSynapsesToSubsetOfWinnerCells")
  (let* ( (tm (make-tm 32 1
              `[basal-input-size . ,(* 32 1)]
              `[activation-threshold . 3]
              `[min-threshold . 1]
              `[sample-size . 3]
              `[predicted-segment-decrement . 0]))
        (matching-segment (create-basal-segment tm 4)))
  (create-basal-synapses tm matching-segment '(0))
  (seg-synapses-set! matching-segment (vector (tm:synapse 0 5000)))
  (compute tm '(0 1 2 3) #t)
  [expect ([tm:get-winner-cells tm] '(0 1 2 3))]
  (compute tm '(4) #t)
  (let* ( (synapses (seg-synapses matching-segment)))
    [expect ((vector-length synapses) 3)]
    (when (= (vector-length synapses) 3)
      [expect ((syn-perm (vector-ref synapses 1)) 2100)
              ((not (member (syn-prex (vector-ref synapses 1)) '(1 2 3))) #f)
              ((syn-perm (vector-ref synapses 2)) 2100)
              ((not (member (syn-prex (vector-ref synapses 2)) '(1 2 3))) #f)])))
                                                                                          ;
(test "MatchingSegmentAddSynapsesToAllWinnerCells")
  (let* ( (tm (make-tm 32 1
              `[basal-input-size . ,(* 32 1)]
              `[activation-threshold . 3]
              `[min-threshold . 1]
              `[sample-size . 3]
              `[predicted-segment-decrement . 0]))
        (matching-segment (create-basal-segment tm 4)))
  (create-basal-synapses tm matching-segment '(0))
  (seg-synapses-set! matching-segment (vector (tm:synapse 0 5000)))
  (compute tm '(0 1) #t)
  [expect ([tm:get-winner-cells tm] '(0 1))]
  (compute tm '(4) #t)
  (let ((synapses (seg-synapses matching-segment)))
    [expect ((vector-length synapses) 2)]
    (when (= (vector-length synapses) 2)
      [expect ((syn-perm (vector-ref synapses 1)) 2100)
              ((not (member (syn-prex (vector-ref synapses 1)) '(1 2 3))) #f)])))
                                                                                          ;
(test "ActiveSegmentGrowSynapsesAccordingToPotentialOverlap")
  (let* ( (tm (make-tm 32 1
              `[basal-input-size . ,(* 32 1)]
              `[activation-threshold . 2]
              `[min-threshold . 1]
              `[sample-size . 4]
              `[predicted-segment-decrement . 0]))
        (active-segment (create-basal-segment tm 5)))
  (create-basal-synapses tm active-segment '(0 1 2))
  (seg-synapses-set! active-segment (vector (tm:synapse 0 5000) (tm:synapse 1 5000) (tm:synapse 2 2000)))
  (compute tm '(0 1 2 3 4) #t)
  [expect ([tm:get-winner-cells tm] '(0 1 2 3 4))]
  (compute tm '(5) #t)
  (let* ( (synapses (seg-synapses active-segment)))
    [expect ((vector-length synapses) 4)
            ((syn-perm (vector-ref synapses 3)) 2100)
            ((not (member (syn-prex (vector-ref synapses 3)) '(3 4))) #f)]))
                                                                                          ;
  ;; flush any test failures
  (when failures
    (display tests)
    (newline))
  (flush-output-port (current-output-port))


)