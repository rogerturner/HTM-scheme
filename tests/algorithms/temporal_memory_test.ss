#!r6rs

;; === HTM-scheme Temporal Memory Test Copyright 2017 Roger Turner. ===
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Based on code from Numenta Platform for Intelligent Computing (NuPIC) ;;
  ;; which is Copyright (C) 2014-2016, Numenta, Inc.                       ;;
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

  ;; Translated from numenta nupic.../temporal_memory_test.py, specialized for
  ;; testing apical_tiebreak_sequence_memory --
  ;; see comments there for descriptions of functions and parameters.
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.
                                                                                            ;
(import 
  (rnrs)
  (HTM-scheme HTM-scheme algorithms htm_prelude)
  (HTM-scheme HTM-scheme algorithms htm_concept)
  (prefix (HTM-scheme HTM-scheme algorithms apical_tiebreak_temporal_memory) attm:)
  (prefix (HTM-scheme HTM-scheme algorithms apical_tiebreak_sequence_memory) atsm:))
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
(define (make-test-tm . options)
  (random-seed! 42)
  (atsm:make-tm (append 
    (if (pair? options) (car options) '())
    `(
      [column-count         . 32]
      [cells-per-column     . 4]
      [activation-threshold . 3]
      [min-threshold        . 2]
      [sample-size          . 3]
      [basal-predicted-segment-decrement  . 0]
      [apical-predicted-segment-decrement . 0]
      [initial-permanence . ,(perm 0.5)] ;; so that grow-synapses creates connected synapses
    ))))
                                                                                            ;
(define (create-basal-segment tm segx)
  (attm:test:create-segment tm segx (attm:test:tm-basal-connections tm)))
                                                                                            ;
(define (create-basal-synapses tm seg pre-cells)
  (attm:test:grow-synapses tm seg (length pre-cells) pre-cells (attm:test:tm-basal-pre-index tm)))
                                                                                            ;
(define-syntax synapses: (identifier-syntax fxvector))
                                                                                            ;
(define (synapses->vector syns)          ;; Synapses -> (vectorof Synapse)
  (let ((vec (make-vector (synapses-length syns))))
    (do ((vx 0 (add1 vx))) ((= vx (synapses-length syns)) vec)
      (vector-set! vec vx (synapses-ref syns vx)))))

(test "ActivateCorrectlyPredictiveCells")
  (let* ( (tm (make-test-tm))
        (active-segment (create-basal-segment tm 4)))
  (create-basal-synapses tm active-segment '(0 1 2 3))
  (atsm:compute tm '(0) #t)
  (atsm:compute tm '(1) #t)
  [expect ([atsm:get-predicted-cells tm] '(4))
          ([atsm:get-active-cells    tm] '(4))])
                                                                                            ;
(test "BurstUnpredictedColumns")
  (let* ( (tm (make-test-tm)))
  (atsm:compute tm '(0) #t)
  [expect ([atsm:get-predicted-cells tm] '())
          ([atsm:get-active-cells    tm] '(0 1 2 3))])
                                                                                            ;
(test "ZeroActiveColumns")
  (let* ( (tm (make-test-tm `(
                [predicted-segment-decrement . ,(perm 0.02)])))
          (segment (create-basal-segment tm 4)))
  (create-basal-synapses tm segment '(0 1 2 3))
  (atsm:compute tm '(0) #t)
  [expect ((null? [atsm:get-active-cells tm]) #f)
          ((null? [attm:get-winner-cells tm])  #f)]
  (atsm:compute tm '() #t)
  [expect ((null? [atsm:get-active-cells tm]) #t)
          ((null? [attm:get-winner-cells tm])  #t)])
                                                                                            ;
(test "PredictedActiveCellsAreAlwaysWinners")
  (let* ( (tm (make-test-tm))
        (active-segment-1 (create-basal-segment tm 4))
        (active-segment-2 (create-basal-segment tm 6)))
  (create-basal-synapses tm active-segment-1 '(0 1 2))
  (create-basal-synapses tm active-segment-2 '(0 1 2))
  (atsm:compute tm '(0) #f)
  (atsm:compute tm '(1) #f)
  [expect ([attm:get-winner-cells tm]  '(4 6))])
                                                                                            ;
(test "ChooseOneWinnerCellInBurstingColumn")
  (let* ( (tm (make-test-tm)))
  (atsm:compute tm '(0) #f)
  [expect ((length [attm:get-winner-cells tm]) 1)]
  (unless (null? [attm:get-winner-cells tm])
    [expect ((not (member (car [attm:get-winner-cells tm]) '(0 1 2 3))) #f)]))
                                                                                            ;
(test "ReinforceCorrectlyActiveSegments")
  (let* ( (tm (make-test-tm `(
                [sample-size                 . 4]
                [permanence-decrement        . ,(perm 0.08)]
                [basal-predicted-segment-decrement . ,(perm 0.02)])))
          (active-segment (create-basal-segment tm 5)))
  (create-basal-synapses tm active-segment '(0 1 2 81))
  (atsm:compute tm '(0) #t)
  (atsm:compute tm '(1) #t)
  [expect ( (vector-map syn-perm (synapses->vector [seg-synapses active-segment]))
            (vector (perm .6) (perm .6) (perm .6) (perm .42)))])
                                                                                            ;
(test "ReinforceSelectedMatchingSegmentInBurstingColumn")
  (let* ( (tm (make-test-tm `(
                [permanence-decrement . ,(perm 0.08)]
                [initial-permanence   . ,(perm 0.3)])))
          (selected-matching-segment (create-basal-segment tm 4))
          (other-matching-segment    (create-basal-segment tm 5)))
  (create-basal-synapses tm selected-matching-segment '(0 1 2 81))
  (create-basal-synapses tm other-matching-segment    '(0 1 81))
  (atsm:compute tm '(0) #t)
  (atsm:compute tm '(1) #t)
  [expect ( (vector-map syn-perm (synapses->vector [seg-synapses selected-matching-segment]))
            (vector (perm .4) (perm .4) (perm .4) (perm .22)))])
                                                                                          ;
(test "NoChangeToNonselectedMatchingSegmentsInBurstingColumn")
  (let* ( (tm (make-test-tm `(
                [permanence-decrement . ,(perm 0.08)]
                [initial-permanence   . ,(perm 0.3)])))
          (selected-matching-segment (create-basal-segment tm 4))
          (other-matching-segment    (create-basal-segment tm 5)))
  (create-basal-synapses tm selected-matching-segment '(0 1 2 81))
  (create-basal-synapses tm other-matching-segment    '(0 1 81))
  (atsm:compute tm '(0) #t)
  (atsm:compute tm '(1) #t)
  [expect ( (vector-map syn-perm (synapses->vector [seg-synapses other-matching-segment]))
            (vector (perm .3) (perm .3) (perm .3)))])
                                                                                          ;
(test "NoChangeToMatchingSegmentsInPredictedActiveColumn")
  (let* ( (tm (make-test-tm))
          (active-segment (create-basal-segment tm 4))
          (matching-segment-on-same-cell  (create-basal-segment tm 4))
          (matching-segment-on-other-cell (create-basal-segment tm 5))
          (synapses (synapses: (make-syn 0 (perm .3)) (make-syn 1 (perm .3)))))
  (create-basal-synapses tm active-segment '(0 1 2 3))
  (create-basal-synapses tm matching-segment-on-same-cell '(0 1))
  (seg-synapses-set! matching-segment-on-same-cell synapses)
  (create-basal-synapses tm matching-segment-on-other-cell '(0 1))
  (seg-synapses-set! matching-segment-on-other-cell synapses)
  (atsm:compute tm '(0) #t)
  (atsm:compute tm '(1) #t)
  [expect ( [atsm:get-predicted-cells tm] '(4))
          ( (vector-map syn-perm (synapses->vector [seg-synapses matching-segment-on-same-cell]))
            (vector (perm .3) (perm .3)))
          ( (vector-map syn-perm (synapses->vector [seg-synapses matching-segment-on-other-cell]))
            (vector (perm .3) (perm .3)))])
                                                                                          ;
(test "NoNewSegmentIfNotEnoughWinnerCells")
  (let* ( (tm (make-test-tm `(
                [sample-size                 . 2]))))
  (atsm:compute tm '() #t)
  (atsm:compute tm '(0) #t)
  [expect ([attm:test:tm-next-flatx tm] 0)])
                                                                                          ;
(test "NewSegmentAddSynapsesToSubsetOfWinnerCells")
  (let* ( (tm (make-test-tm `(
                [sample-size        . 2]
                [initial-permanence . ,(perm 0.21)]))))
  (atsm:compute tm '(0 1 2) #t)
  (let ((prev-winner-cells [attm:get-winner-cells tm]))
    [expect ((length prev-winner-cells) 3)]
    (atsm:compute tm '(4) #t)
    [expect ((length [attm:get-winner-cells tm]) 1)]
    (let* ( (segments (vector-ref [attm:test:tm-basal-connections tm] (car [attm:get-winner-cells tm])))
            (synapses (seg-synapses (car segments))))
      [expect ((length segments) 1)
              ((synapses-length synapses) 2)]
      (vector-for-each
        (lambda (synapse)
          [expect ((syn-perm synapse) (perm .21))
                  ((not (member (syn-prex synapse) prev-winner-cells)) #f)])
        (synapses->vector synapses)))))
                                                                                          ;
(test "NewSegmentAddSynapsesToAllWinnerCells")
  (let* ( (tm (make-test-tm `(
                [sample-size        . 4]
                [initial-permanence . ,(perm 0.21)]))))
  (atsm:compute tm '(0 1 2) #t)
  (let ((prev-winner-cells [attm:get-winner-cells tm]))
    [expect ((length prev-winner-cells) 3)]
    (atsm:compute tm '(4) #t)
    [expect ((length [attm:get-winner-cells tm]) 1)]
    (let* ( (segments (vector-ref [attm:test:tm-basal-connections tm] (car [attm:get-winner-cells tm])))
            (synapses (seg-synapses (car segments))))
      [expect ((length segments) 1)
              ((synapses-length synapses) 3)]
      (vector-for-each
        (lambda (synapse)
          [expect ((syn-perm synapse) (perm .21))
                  ((not (member (syn-prex synapse) prev-winner-cells)) #f)])
        (synapses->vector synapses))
      [expect ((vector->list (vector-map syn-prex (synapses->vector synapses))) prev-winner-cells)])))
                                                                                          ;
(test "MatchingSegmentAddSynapsesToSubsetOfWinnerCells")
  (let* ( (tm (make-test-tm `(
                [cells-per-column . 1]
                [basal-input-size . ,(* 32 1)]
                [min-threshold    . 1]
                [initial-permanence . ,(perm 0.21)])))
        (matching-segment (create-basal-segment tm 4)))
  (create-basal-synapses tm matching-segment '(0))
  (seg-synapses-set! matching-segment (synapses: (make-syn 0 (perm .5))))
  (atsm:compute tm '(0 1 2 3) #t)
  [expect ([attm:get-winner-cells tm] '(0 1 2 3))]
  (atsm:compute tm '(4) #t)
  (let* ( (synapses (seg-synapses matching-segment)))
    [expect ((synapses-length synapses) 3)]
    (when (= (synapses-length synapses) 3)
      [expect ((syn-perm (synapses-ref synapses 1)) (perm .21))
              ((not (member (syn-prex (synapses-ref synapses 1)) '(1 2 3))) #f)
              ((syn-perm (synapses-ref synapses 2)) (perm .21))
              ((not (member (syn-prex (synapses-ref synapses 2)) '(1 2 3))) #f)])))
                                                                                          ;
(test "MatchingSegmentAddSynapsesToAllWinnerCells")
  (let* ( (tm (make-test-tm `(
                [cells-per-column . 1]
                [basal-input-size . ,(* 32 1)]
                [min-threshold    . 1]
                [initial-permanence . ,(perm 0.21)])))
          (matching-segment (create-basal-segment tm 4)))
  (create-basal-synapses tm matching-segment '(0))
  (seg-synapses-set! matching-segment (synapses: (make-syn 0 (perm .5))))
  (atsm:compute tm '(0 1) #t)
  [expect ([attm:get-winner-cells tm] '(0 1))]
  (atsm:compute tm '(4) #t)
  (let ((synapses (seg-synapses matching-segment)))
    [expect ((synapses-length synapses) 2)]
    (when (= (synapses-length synapses) 2)
      [expect ((syn-perm (synapses-ref synapses 1)) (perm .21))
              ((not (member (syn-prex (synapses-ref synapses 1)) '(1 2 3))) #f)])))
                                                                                          ;
(test "ActiveSegmentGrowSynapsesAccordingToPotentialOverlap")
  (let* ( (tm (make-test-tm `(
                [cells-per-column     . 1]
                [basal-input-size . ,(* 32 1)]
                [activation-threshold . 2]
                [min-threshold . 1]
                [sample-size . 4]
                [initial-permanence . ,(perm 0.21)])))
          (active-segment (create-basal-segment tm 5)))
  (create-basal-synapses tm active-segment '(0 1 2))
  (seg-synapses-set! active-segment
    (synapses: (make-syn 0 (perm .5)) (make-syn 1 (perm .5)) (make-syn 2 (perm .2))))
  (atsm:compute tm '(0 1 2 3 4) #t)
  [expect ([attm:get-winner-cells tm] '(0 1 2 3 4))]
  (atsm:compute tm '(5) #t)
  (let* ( (synapses (seg-synapses active-segment)))
    [expect ((synapses-length synapses) 4)
            ((syn-perm (synapses-ref synapses 3)) (perm .21))
            ((not (member (syn-prex (synapses-ref synapses 3)) '(3 4))) #f)]))
                                                                                          ;
  ;; flush any test failures
  (when #t ;failures
    (display tests)
    (newline))
  (flush-output-port (current-output-port))
