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
        (only (chezscheme) time)
        (rnrs)
        (HTM-scheme HTM-scheme algorithms htm_prelude)
        (HTM-scheme HTM-scheme algorithms apical_tiebreak_sequence_memory))

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
              (break)))))
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
                  (atsm:compute tm (car pattern) #f)
                  (when (positive? i)
                    (assert-equal (set (atsm:get-predicted-cells tm))
                                  (set (atsm:get-active-cells tm))))
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
        `([initial-permanence   . ,(tm:permanence 0.2)]
          [connected-permanence . ,(tm:permanence 0.7)]
          [permanence-increment . ,(tm:permanence 0.2)]))
(test-B "B8: like B7 but with 32 cellsPerColumn"
        1 N 4 32
        `([initial-permanence   . ,(tm:permanence 0.2)]
          [connected-permanence . ,(tm:permanence 0.7)]
          [permanence-increment . ,(tm:permanence 0.2)]))
    )))))
  
;; --- Helper functions ---
                                                                                            ;
(define (init cpc . options)
  (apply atsm:construct 
    (append (list n cpc)
            (car options)
          `([initial-permanence          . ,(tm:permanence 0.5)]
            [connected-permanence        . ,(tm:permanence 0.6)]
            [min-threshold               . 25]
            [sample-size                 . 30]
            [permanence-increment        . ,(tm:permanence 0.1)]
            [permanence-decrement        . ,(tm:permanence 0.02)]
            [predicted-segment-decrement . ,(tm:permanence 0.08)] ;;ie basal-p-s-d
            [activation-threshold        . 25]
            [seed                        . 42]
            [max-synapses-per-segment    . -1]
            ))))
                                                                                            ;
(define (random-pattern)
  ;; sorted random selection of w integers in range 0..n-1
  (list-sort < (vector->list (vector-sample (build-vector n id) w))))
                                                                                            ;
(define (set l)
  (unique = l))
  