;; === HTM-scheme Apical Tiebreak Test Base Copyright 2017 Roger Turner. ===
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

  ;; Translated from numenta htmresearch/.../apical_tiebreak_test_base.py --
  ;; see comments there for descriptions of functions and parameters.
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.

#!chezscheme 
                                                                                            ;
(library-directories "../../../HTM-scheme/algorithms/")
                                                                                            ;
(import (except (chezscheme) add1 make-list random)
        (htm_prelude)
        (apical_tiebreak_pair_memory))

(define apical-input-size 1000)
(define basal-input-size  1000)
(define column-count      2048)
(define w                 40)
(define cells-per-column  32)

;; --- Helper functions ---
                                                                                            ;
(define (init)                           ;; -> TM
  (atpm:construct column-count cells-per-column `(
    [basal-input-size            . ,basal-input-size]
    [apical-input-size           . ,apical-input-size]
    [initial-permanence          . ,(tm:permanence 0.5)]
    [connected-permanence        . ,(tm:permanence 0.6)]
    [min-threshold               . 25]
    [sample-size                 . 30]
    [permanence-increment        . ,(tm:permanence 0.1)]
    [permanence-decrement        . ,(tm:permanence 0.02)]
    [basal-predicted-segment-decrement . ,(tm:permanence 0.0)]
    [activation-threshold        . 25]
    [seed                        . 42]
    )))
                                                                                            ;
(define (cellx->colx cellx)              ;; CellX -> ColX
  (quotient cellx cells-per-column))
                                                                                            ;
(define (list-or l1 l2)                  ;; {Nat} {Nat} -> {Nat}
  (bitwise->list (bitwise-ior (list->bitwise l1) (list->bitwise l2))))
                                                                                            ;
(define (list-difference l1 l2)          ;; {Nat} {Nat} -> {Nat}
  (bitwise->list (bitwise-and (list->bitwise l1) (bitwise-not (list->bitwise l2)))))
                                                                                            ;
(define (get-bursting-columns tm)        ;; TM -> (listof ColX)
  (let ((predicted (unique = (map cellx->colx (atpm:get-predicted-cells tm))))
        (active    (unique = (map cellx->colx (atpm:get-active-cells tm)))))
    (list-difference active predicted)))
                                                                                            ;
(define (random-pattern size)            ;; Nat -> {Nat}
  ;; sorted random selection of w integers in range 0..size-1
  (list-sort < (vector->list (vector-sample (build-vector size id) w))))

;; --- Unit testing ---
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

(test "BasalInputRequiredForPredictions")
  ;; Learn A for basalInput1, apicalInput1.
  ;; Now observe A with apicalInput1 but no basal input. It should burst.
  (let ((tm  (init))
        (active-columns (random-pattern column-count))
        (basal-input    (random-pattern basal-input-size))
        (apical-input   (random-pattern apical-input-size)))
    (do ((x 0 (add1 x))) ((= x 3))
      (atpm:compute tm active-columns basal-input apical-input '() '() #t))
    (atpm:compute tm active-columns '() apical-input '() '() #f)
    (expect ((get-bursting-columns tm) active-columns)))
                                                                                            ;
(test "BasalPredictionsWithoutApical")
  ;; Learn A for two contexts:
  ;;  - basalInput1, apicalInput1
  ;;  - basalInput2, apicalInput2
  ;; Now observe A with a union of basalInput1 and basalInput2, and no apical
  ;; input. It should predict both contexts.
  (let ((tm  (init))
        (active-columns  (random-pattern column-count))
        (basal-input-1   (random-pattern basal-input-size))
        (basal-input-2   (random-pattern basal-input-size))
        (apical-input-1  (random-pattern apical-input-size))
        (apical-input-2  (random-pattern apical-input-size)))
    (do ((x 0 (add1 x))) ((= x 3))
      (atpm:compute tm active-columns basal-input-1 apical-input-1 '() '() #t)
      (set! active-cells-1 (atpm:get-active-cells tm))
      (atpm:compute tm active-columns basal-input-2 apical-input-2 '() '() #t)
      (set! active-cells-2 (atpm:get-active-cells tm)))
    (atpm:compute tm active-columns
                (list-or basal-input-1 basal-input-2)
                '() '() '() #f)
    (expect ((atpm:get-active-cells tm) (list-or active-cells-1 active-cells-2))))
                                                                                            ;
(test "ApicalNarrowsThePredictions")
  ;; Learn A for two contexts:
  ;;  - basalInput1, apicalInput1
  ;;  - basalInput2, apicalInput2
  ;; Now observe A with a union of basalInput1 and basalInput2, and apicalInput1.
  ;; It should only predict one context.
  (let ((tm  (init))
        (active-columns  (random-pattern column-count))
        (basal-input-1   (random-pattern basal-input-size))
        (basal-input-2   (random-pattern basal-input-size))
        (apical-input-1  (random-pattern apical-input-size))
        (apical-input-2  (random-pattern apical-input-size))
        (active-cells-1  #f)
        (active-cells-2  #f))
    (do ((x 0 (add1 x))) ((= x 3))
      (atpm:compute tm active-columns basal-input-1 apical-input-1 '() '() #t)
      (set! active-cells-1 (atpm:get-active-cells tm))
      (atpm:compute tm active-columns basal-input-2 apical-input-2 '() '() #t)
      (set! active-cells-2 (atpm:get-active-cells tm)))
    (atpm:compute tm active-columns
                (list-or basal-input-1 basal-input-2)
                apical-input-1 '() '() #f)
    (expect ((atpm:get-active-cells tm) active-cells-1)))
                                                                                            ;
(test "UnionOfFeedback")
  ;; Learn A for three contexts:
  ;;  - basalInput1, apicalInput1
  ;;  - basalInput2, apicalInput2
  ;;  - basalInput3, apicalInput3
  ;; Now observe A with a union of all 3 basal inputs, and a union of
  ;; apicalInput1 and apicalInput2. It should predict 2 of the 3 contexts.
  (let ((tm  (init))
        (active-columns  (random-pattern column-count))
        (basal-input-1   (random-pattern basal-input-size))
        (basal-input-2   (random-pattern basal-input-size))
        (basal-input-3   (random-pattern basal-input-size))
        (apical-input-1  (random-pattern apical-input-size))
        (apical-input-2  (random-pattern apical-input-size))
        (apical-input-3  (random-pattern apical-input-size))
        (active-cells-1  #f)
        (active-cells-2  #f)
        (active-cells-3  #f))
    (do ((x 0 (add1 x))) ((= x 3))
      (atpm:compute tm active-columns basal-input-1 apical-input-1 '() '() #t)
      (set! active-cells-1 (atpm:get-active-cells tm))
      (atpm:compute tm active-columns basal-input-2 apical-input-2 '() '() #t)
      (set! active-cells-2 (atpm:get-active-cells tm))
      (atpm:compute tm active-columns basal-input-3 apical-input-3 '() '() #t)
      (set! active-cells-3 (atpm:get-active-cells tm)))
    (atpm:compute tm active-columns
                (list-or (list-or basal-input-1 basal-input-2) basal-input-3)
                (list-or apical-input-1 apical-input-2) '() '() #f)
    (expect ((atpm:get-active-cells tm) (list-or active-cells-1 active-cells-2))))
    
(display tests) (newline)
