#!r6rs

;; === HTM-scheme Sequence Memory Test Base Copyright 2018 Roger Turner. ===
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

  ;; Translated from numenta htmresearch/.../sequence_memory_test_base.py --
  ;; see comments there for descriptions of functions and parameters;
  ;; specialized to apical_tiebreak_sequence_memory.
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.

(import
  (rnrs)
  (HTM-scheme HTM-scheme algorithms htm_prelude)
  (HTM-scheme HTM-scheme algorithms htm_concept)
  (prefix (HTM-scheme HTM-scheme algorithms column_pooler) atsm:))

(define n            2048)
(define w              40)
(define feedback-size 400)
(define N             100)

(define run 
  (lambda () (call/cc (lambda (break)  ;; exit to repl on failure
  (let* 
    ( (list-summary 
        (lambda (l)
          (append (list (length l) ': ) (take 6 l) (list '- ) 
                  (reverse (take 6 (reverse l))))))
      (show
        (lambda (x)
          (display x)
          (flush-output-port (current-output-port))))
      (assert-equal 
        (lambda (x y)
          (unless (equal? x y)
            (begin
              (for-each display `(
                "assertion failed\n"
                ,(list-summary x) #\newline
                ,(list-summary y) #\newline))
              #;(break)))))
#|
  Input Sequence: We train with M input sequences, each consisting of N random
  patterns. Each pattern consists of a random number of bits on. The number of
  1's in each pattern should be between 38 and 40 columns.

  Each input pattern can optionally have an amount of spatial noise represented
  by X, where X is the probability of switching an on bit with a random bit.

  Training: The ETM is trained with P passes of the M sequences. There
  should be a reset between sequences. The total number of iterations during
  training is P*N*M.

  Testing: Run inference through the same set of sequences, with a reset before
  each sequence. For each sequence the system should accurately predict the
  pattern at the next time step up to and including the N-1'st pattern. The
  number of predicted inactive cells at each time step should be reasonably low.

  We can also calculate the number of synapses that should be
  learned. We raise an error if too many or too few were learned.
|#
      (test-B 
        (lambda (title M N P cpc options)
          (show title)
          (let ((tm        (init cpc options))
                (sequences (build-list M (lambda (_)
                             (build-list N (lambda (_) (random-pattern)))))))
            ;; Learn
            (do ((_ 0 (add1 _))) ((= _ P))
              (for-each
                (lambda (sequence)
                  (for-each
                    (lambda (pattern)
                      (atsm:compute tm pattern #t))
                    sequence)
                  (atsm:reset tm)
                  (show " ."))
                sequences))
            ;; Predict
            (for-each
              (lambda (sequence)
                (let loop ((i 0) (pattern sequence))
                  (let ((predicted-cells (set (atsm:get-active-cells tm))))
                    (atsm:compute tm (car pattern) #f)
                    (when (positive? i)
                      (assert-equal predicted-cells ;(set (atsm:get-predicted-cells tm))
                                    (set (atsm:get-active-cells tm)))))
                  (when (< i (- N 1))
                    (loop (add1 i) (cdr pattern))))
                (atsm:reset tm)
                (show " ."))
              sequences)
            (show " ok\n")))))
          
(test-B "B1: Basic sequence learner.  M=1, N=100, P=2"
        1 N 2 1 '())
(test-B "B3: M=1, N=300, P=2. (See how high we can go with N)"
        1 (* 3 N) 2 1 '())
(test-B "B4: M=3, N=300, P=2. (See how high we can go with N*M)"
        3 (* 3 N) 2 1 '())
(test-B "B5: like B1 but with 32 cellsPerColumn"
        1 N 2 32 '())
(test-B "B6: like B4 but with 32 cellsPerColumn"
        3 (* 3 N) 2 32 '())
(test-B "B7: like B1 but with slower learning"
        1 N 4 1 
        `([initial-permanence   . ,(perm 0.2)]
          [connected-permanence . ,(perm 0.7)]
          [permanence-increment . ,(perm 0.2)]))
(test-B "B8: like B7 but with 32 cellsPerColumn"
        1 N 4 32
        `([initial-permanence   . ,(perm 0.2)]
          [connected-permanence . ,(perm 0.7)]
          [permanence-increment . ,(perm 0.2)]))
    )))))
  
;; --- Helper functions ---
                                                                                            ;
(define (init cpc options)
  (atsm:make-cp
    (append 
      options `(
      [cell-count                  . ,n]
      [cells-per-column            . ,cpc]
      [initial-permanence          . ,(perm 0.5)]
      [connected-permanence        . ,(perm 0.6)]
      [min-threshold               . 25]
      [sample-size                 . 30]
      [permanence-increment        . ,(perm 0.1)]
      [permanence-decrement        . ,(perm 0.02)]
      [basal-predicted-segment-decrement . ,(perm 0.08)]
      [activation-threshold        . 25]
      [seed                        . 42]
      [max-synapses-per-segment    . -1]))))
                                                                                            ;
(define (random-pattern)
  ;; sorted random selection of w integers in range 0..n-1
  (list-sort < (vector->list (vector-sample (build-vector n id) w))))
                                                                                            ;
(define (set l)
  (unique = l))
  
(run)
  