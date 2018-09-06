#!r6rs

;; === HTM-scheme Apical Tiebreak Sequences Test Base Copyright 2017 Roger Turner. ===
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

  ;; Translated from numenta htmresearch/.../apical_tiebreak_sequences_test_base.py --
  ;; see comments there for descriptions of functions and parameters.
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.

(import
        ;(only (racket) time-apply)
        (only (chezscheme) time)
        (rnrs)
        (HTM-scheme HTM-scheme algorithms htm_prelude)
        (HTM-scheme HTM-scheme algorithms apical_tiebreak_sequence_memory))

(define column-count      2048)
(define w                   40)
(define apical-input-size 1000)
(define cells-per-column    32)

;; --- Helper functions ---
                                                                                            ;
(define (init . options)                           ;; -> TM
  (apply atsm:construct 
    (append (list column-count cells-per-column)
            options
          `([apical-input-size           . ,apical-input-size]
            [initial-permanence          . ,(tm:permanence 0.5)]
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
(define (cellx->colx cellx)              ;; CellX -> ColX
  (div cellx cells-per-column))
                                                                                            ;
(define (list-or l1 l2)                  ;; {Nat} {Nat} -> {Nat}
  (bitwise->list (bitwise-ior (list->bitwise l1) (list->bitwise l2))))
                                                                                            ;
(define (list-difference l1 l2)          ;; {Nat} {Nat} -> {Nat}
  (bitwise->list (bitwise-and (list->bitwise l1) (bitwise-not (list->bitwise l2)))))
                                                                                            ;
(define (get-bursting-columns tm)        ;; TM -> (listof ColX)
  (let ((predicted (unique = (map cellx->colx (atsm:get-predicted-cells tm))))
        (active    (unique = (map cellx->colx (atsm:get-active-cells tm)))))
    (list-difference active predicted)))
                                                                                            ;
(define (random-pattern size)            ;; Nat -> {Nat}
  ;; sorted random selection of w integers in range 0..size-1
  (list-sort < (vector->list (vector-sample (build-vector size id) w))))
                                                                                            ;
(define (random-column-pattern)        ;; Nat -> {Nat}
  (random-pattern column-count))
                                                                                            ;
(define (random-apical-pattern)        ;; Nat -> {Nat}
  (random-pattern apical-input-size))
  
(define (filter-cells-by-column cells columns)
  (filter
    (lambda (cell)
      (member (div cell cells-per-column) columns))
    cells))
    
(define (set l)   ;; {X} -> {X}
  (unique = l))
  
;; --- Unit testing ---
                                                                                            ;
(define tests "")
(define failures #f)
                                                                                            ;
(define (test name)
  (set! tests (string-append tests "\n" name "\n")))

(define (dl l)
  (append (list (length l) ': )(take 6 l) (list '- ) (reverse (take 6 (reverse l)))))
                                                                                            ;
(define-syntax expect                    ;; ((X ... -> Y) X ...) Y -> [error]
  ;; check that function application(s) to arguments match expected values
  (lambda (x)                            
    (syntax-case x ()                                  ;; [expect ([fn args] expected ) ... ]
      [ (_ (expr expected) ...)                        ;; expr matches [fn args]
        #'(begin (let ((result expr))                  ;; eval expr just once
                   (if (equal? result expected)
                      (begin
                        (for-each display `(,tests expr " ok" #\newline))
                        (set! tests "")
                        )
                      (begin 
                        (for-each display `(,tests #\newline expr " **" #\newline
                          "  expected: " ,(dl expected) #\newline 
                          "  returned: " ,(dl result) #\newline))
                        (set! failures #t)
                        (set! tests "")))) ...)])))
                                                                                            ;
(define e-cells #f) 
(define y-cells #f)
                                                                                            ;

(define run (lambda ()

(test "SequenceMemory_BasalInputRequiredForPredictions")
  ;; Learn ABCDE with F1.
  ;; Reset, then observe B with F1.
  ;; It should burst, despite the fact that the B cells have apical support.
  (let* ((tm       (init))
        (r        random-column-pattern)
        (abcde    (cons (r) (cons (r) (cons (r) (cons (r) (cons (r) '()))))))
        (feedback (random-apical-pattern)))
    (do ((_ 0 (add1 _))) ((= _ 4))
      (atsm:reset tm)
      (for-each
        (lambda (pattern)
          (atsm:compute tm pattern feedback #t))
        abcde))
    (atsm:reset tm)
    (atsm:compute tm (list-ref abcde 1) feedback #f)
    (expect ( (atsm:get-predicted-cells tm)  '() )
            ( (get-bursting-columns tm) (list-ref abcde 1) )))
                                                                                            ;
(test "SequenceMemory_BasalPredictionsWithoutFeedback")
  ;; Train on ABCDE with F1, XBCDY with F2.
  ;; Test with BCDE. Without feedback, two patterns are predicted.
  (let* ((tm       (init))
        (r         random-column-pattern)
        (bcd       (cons (r) (cons (r) (cons (r) '()))))
        (abcde     (cons (random-column-pattern) 
                     (append bcd (list (random-column-pattern)))))
        (xbcdy     (cons (random-column-pattern) 
                     (append bcd (list (random-column-pattern)))))
        (feedback1 (random-apical-pattern))
        (feedback2 (random-apical-pattern)))
    ;; First learn the sequences without feedback. We need to let it work through
    ;; the common subsequence, choosing new cell SDRs for elements later in the
    ;; sequence, before allowing it to grow apical segments.
    (do ((_ 0 (add1 _))) ((= _ 10))
      (atsm:reset tm)
      (for-each
        (lambda (pattern)
          (atsm:compute tm pattern '() #t))
        xbcdy)
      (atsm:reset tm)
      (for-each
        (lambda (pattern)
          (atsm:compute tm pattern '() #t))
        abcde))
    ;; Learn the apical connections
    (do ((_ 0 (add1 _))) ((= _ 2))
      (atsm:reset tm)
      (for-each
        (lambda (pattern)
          (atsm:compute tm pattern feedback2 #t))
        xbcdy)
      (set! y-cells (set (atsm:get-active-cells tm)))
      (atsm:reset tm)
      (for-each
        (lambda (pattern)
          (atsm:compute tm pattern feedback1 #t))
        abcde)
      (set! e-cells (set (atsm:get-active-cells tm))))
    ;; Test    
    (atsm:reset tm)
    (for-each
      (lambda (pattern)
        (atsm:compute tm pattern '() #f))
      (cdr abcde))
    ;; The E cells should be active, and so should any Y cells that happen to be
    ;; in a minicolumn shared between E and Y.
    
    (let ((expected-active
            (list-or e-cells (set (filter-cells-by-column y-cells (list-ref abcde 4))))))
      (expect ((set (atsm:get-active-cells tm)) expected-active )
              ((set (atsm:get-predicted-cells tm)) (list-or e-cells y-cells)))))

(test "SequenceMemory_FeedbackNarrowsThePredictions")
  ;; Train on ABCDE with F1, XBCDY with F2.
  ;; Test with BCDE with F1. One pattern is predicted.
  (let* ( (tm        (init))
          (r         random-column-pattern)
          (bcd       (cons (r) (cons (r) (cons (r) '()))))
          (abcde     (cons (random-column-pattern) 
                       (append bcd (list (random-column-pattern)))))
          (xbcdy     (cons (random-column-pattern) 
                       (append bcd (list (random-column-pattern)))))
          (feedback1 (random-apical-pattern))
          (feedback2 (random-apical-pattern)))
  ;; First learn the sequences without feedback. We need to let it work through
  ;; the common subsequence, choosing new cell SDRs for elements later in the
  ;; sequence, before allowing it to grow apical segments.
    (do ((_ 0 (add1 _))) ((= _ 10))
      (atsm:reset tm)
      (for-each
        (lambda (pattern)
          (atsm:compute tm pattern '() #t))
        abcde)
      (atsm:reset tm)
      (for-each
        (lambda (pattern)
          (atsm:compute tm pattern '() #t))
        xbcdy))
    ;; Learn the apical connections
    (do ((_ 0 (add1 _))) ((= _ 2))
      (atsm:reset tm)
      (for-each
        (lambda (pattern)
          (atsm:compute tm pattern feedback1 #t))
        abcde)
      (set! e-cells (atsm:get-active-cells tm))
      (atsm:reset tm)
      (for-each
        (lambda (pattern)
          (atsm:compute tm pattern feedback2 #t))
        xbcdy)
      (set! y-cells (atsm:get-active-cells tm)))
    ;; Test    
    (atsm:reset tm)
    (for-each
      (lambda (pattern)
        (atsm:compute tm pattern feedback1 #f))
      (cdr abcde))
    (expect ((atsm:get-active-cells tm) e-cells )
            ((atsm:get-predicted-cells tm) e-cells)))

(test "SequenceMemory_IncorrectFeedbackLeadsToBursting")
  ;; Train on ABCDE with F1, XBCDY with F2.
  ;; Test with BCDE with F2. E should burst.
  (let* ((tm        (init))
        (r         random-column-pattern)
        (bcd       (cons (r) (cons (r) (cons (r) '()))))
        (abcde     (cons (random-column-pattern) 
                     (append bcd (list (random-column-pattern)))))
        (xbcdy     (cons (random-column-pattern) 
                     (append bcd (list (random-column-pattern)))))
        (feedback1 (random-apical-pattern))
        (feedback2 (random-apical-pattern))
        (e-cells #f) (y-cells #f))
  ;; First learn the sequences without feedback. We need to let it work through
  ;; the common subsequence, choosing new cell SDRs for elements later in the
  ;; sequence, before allowing it to grow apical segments.
    (do ((_ 0 (add1 _))) ((= _ 10))
      (atsm:reset tm)
      (for-each
        (lambda (pattern)
          (atsm:compute tm pattern '() #t))
        abcde)
      (atsm:reset tm)
      (for-each
        (lambda (pattern)
          (atsm:compute tm pattern '() #t))
        xbcdy))
    ;; Learn the apical connections
    (do ((_ 0 (add1 _))) ((= _ 2))
      (atsm:reset tm)
      (for-each
        (lambda (pattern)
          (atsm:compute tm pattern feedback1 #t))
        abcde)
      (set! e-cells (atsm:get-active-cells tm))
      (atsm:reset tm)
      (for-each
        (lambda (pattern)
          (atsm:compute tm pattern feedback2 #t))
        xbcdy)
      (set! y-cells (atsm:get-active-cells tm)))
    ;; Test    
    (atsm:reset tm)
    (for-each
      (lambda (pattern)
        (atsm:compute tm pattern feedback2 #f))
      (cdr abcde))
    (expect ((atsm:get-predicted-cells tm) y-cells )
            ((get-bursting-columns tm) 
               (list-difference (list-ref abcde 4) (list-ref xbcdy 4)))))

(test "SequenceMemory_UnionOfFeedback")
  ;; Train on ABCDE with F1, XBCDY with F2, MBCDN with F3.
  ;; Test with BCDE with F1 | F2. The last step should predict E and Y.
  (let* ((tm        (init))
        (r         random-column-pattern)
        (bcd       (cons (r) (cons (r) (cons (r) '()))))
        (abcde     (cons (random-column-pattern) 
                     (append bcd (list (random-column-pattern)))))
        (xbcdy     (cons (random-column-pattern) 
                     (append bcd (list (random-column-pattern)))))
        (mbcdn     (cons (random-column-pattern) 
                     (append bcd (list (random-column-pattern)))))
        (feedback1 (random-apical-pattern))
        (feedback2 (random-apical-pattern))
        (feedback3 (random-apical-pattern))
        (e-cells #f) (y-cells #f))
  ;; First learn the sequences without feedback. We need to let it work through
  ;; the common subsequence, choosing new cell SDRs for elements later in the
  ;; sequence, before allowing it to grow apical segments.
    (do ((_ 0 (add1 _))) ((= _ 20))
      (atsm:reset tm)
      (for-each
        (lambda (pattern)
          (atsm:compute tm pattern '() #t))
        abcde)
      (atsm:reset tm)
      (for-each
        (lambda (pattern)
          (atsm:compute tm pattern '() #t))
        xbcdy)
      (atsm:reset tm)
      (for-each
        (lambda (pattern)
          (atsm:compute tm pattern '() #t))
        mbcdn))
    ;; Learn the apical connections
    (do ((_ 0 (add1 _))) ((= _ 2))
      (atsm:reset tm)
      (for-each
        (lambda (pattern)
          (atsm:compute tm pattern feedback1 #t))
        abcde)
      (set! e-cells (atsm:get-active-cells tm))
      (atsm:reset tm)
      (for-each
        (lambda (pattern)
          (atsm:compute tm pattern feedback2 #t))
        xbcdy)
      (set! y-cells (atsm:get-active-cells tm))
      (atsm:reset tm)
      (for-each
        (lambda (pattern)
          (atsm:compute tm pattern feedback3 #t))
        mbcdn))
    ;; Test    
    (atsm:reset tm)
    (for-each
      (lambda (pattern)
        (atsm:compute tm pattern (list-or feedback1 feedback2) #f))
      (cdr abcde))
    (let ((expected-active
            (list-or e-cells (filter-cells-by-column y-cells (list-ref abcde 4)))))
      (expect ((atsm:get-active-cells tm) expected-active )
              ((atsm:get-predicted-cells tm) (list-or e-cells y-cells)))))

(display tests) (newline)

              ))

;(time-apply run '())
;(time (run))
(run)
