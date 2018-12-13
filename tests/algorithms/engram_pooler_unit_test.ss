#!r6rs

;; === HTM-scheme Engram Pooler Unit Test Copyright 2018 Roger Turner. ===
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Based on code from Numenta Platform for Intelligent Computing (NuPIC) ;;
  ;; which is Copyright (C) 2016, Numenta, Inc.                            ;;
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

  ;; Translated from numenta htmresearch/.../column_pooler_unit_test.py --
  ;; see comments there for descriptions of functions and parameters.
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.

(import 
          (rnrs)
          (HTM-scheme HTM-scheme algorithms htm_prelude)
          (HTM-scheme HTM-scheme algorithms htm_concept)
  (prefix (HTM-scheme HTM-scheme algorithms engram_pooler) cp:))
                                                                                            ;
(define tests "")
(define successes 0)
(define failures  0)
                                                                                            ;
(define (test name)
  (newline) (display name) (set! tests "")
  #;(set! tests (string-append tests "\n  " name)))
                                                                                            ;
(define-syntax expect                    ;; ((X ... -> Y) X ...) Y -> [error]
  ;; check that function application(s) to arguments match expected values
  (lambda (x)                            
    (syntax-case x ()                                  ;; [expect ([fn args] expected ) ... ]
      [ (_ (expr expected) ...)                        ;; expr matches [fn args]
        #'(begin (let ((result expr))                  ;; eval expr just once
                   (if (equal? result expected)
                      (set! successes (add1 successes))
                      (begin 
                        (for-each display `(,tests #\newline "**" expr #\newline
                          "  expected: " ,expected #\newline 
                          "  returned: " ,result  #\newline))
                        (set! failures (add1 failures))
                        (set! tests "")))) ...)]
      [ (_ expr ...)
        #'(begin (if expr
                    (set! successes (add1 successes))
                    (begin 
                      (for-each display `(,tests #\newline "**" expr " failed\n"))
                      (set! failures (add1 failures))
                      (set! tests ""))) ...)]
                      )))
                                                                                            ;
(define minss 10)
(define maxss 40)
(define (initialize-default-pooler . kwargs)
  (cp:make-cp 
    (append
      kwargs
      `(
        [min-sdr-size . ,minss]
        [max-sdr-size . ,maxss]
        [input-width  . ,(* 2048 8)]
        [cell-count   . 4096]))))
                                                                                            ;
(define (set l)
  (unique = l))
                                                                                            ;
(define (range start limit)              ;; Nat Nat -> {Nat}
  ;; produce list of values from start to (limit-1)
  (build-list (fx- limit start)
    (lambda (i) (+ start i))))
                                                                                            ;
(test "Constructor")                     ;; Create a simple instance and test the constructor
  (let ((pooler (initialize-default-pooler)))
    (expect [(cp:number-of-cells   pooler)  4096]
            [(cp:number-of-inputs  pooler) 16384]
            [(cp:test:number-of-proximal-synapses 
                    pooler (range 0 2048))     0]
            [(cp:test:num-connected-proximal-synapses 
                    pooler (range 0 2048))     0]))
                                                                                            ;
(test "InitialNullInputLearnMode")       ;; Tests with no input in the beginning
  (let ((pooler (initialize-default-pooler)))

    ;; Should be no active cells in beginning
    (expect [(length (cp:get-active-cells pooler))  0])

    ;; After computing with no input should have active cells
    (cp:compute pooler (list) #t)
    (let ((object-sdr1 (set (cp:get-active-cells pooler))))
      (expect [<= minss (length object-sdr1) maxss])
      
      ;; Should be no active cells after reset
      (cp:reset pooler)
      (expect [(length (cp:get-active-cells pooler))  0])
      
      ;; Computing again with no input should lead to different active cells
      (cp:compute pooler (list) #t)
      (let ((object-sdr2 (set (cp:get-active-cells pooler))))
        (expect [<= minss (length object-sdr2) maxss]
                [< (length (intersect1d object-sdr1 object-sdr2)) 5]))))
                                                                                            ;
(test "InitialProximalLearning")         ;; Tests the first few steps of proximal learning
  (let ((pooler (initialize-default-pooler)))
  
    ;; Get initial activity
    (cp:compute pooler (range 0 40) #t)
    (let* ( (object-sdr (set (cp:get-active-cells pooler)))
            (object-sdr-size (length object-sdr)))
      (expect [<= minss object-sdr-size maxss])

      ;; Ensure we've added correct number synapses on the active cells
      ;; Ensure they are all connected
      (expect [(cp:test:num-connected-proximal-synapses pooler (cp:get-active-cells pooler)) (* object-sdr-size 20)])

      ;; As multiple different feedforward inputs come in, the same set of cells should be active
      (cp:compute pooler (range 100 140) #t)
      (expect [(set (cp:get-active-cells pooler)) object-sdr])
      (expect [(cp:test:num-connected-proximal-synapses pooler (cp:get-active-cells pooler)) (* object-sdr-size 40)])

      ;; If there is no feedforward input we should still get the same set of active cells
      (cp:compute pooler (list) #t)
      (expect [(set (cp:get-active-cells pooler)) object-sdr])

      ;; In "learn new object mode", given a familiar feedforward input after reset we should not get the same set of active cells
      (cp:reset pooler)
      (cp:compute pooler (range 0 40) #t)
      (expect [(equal? object-sdr (set (cp:get-active-cells pooler))) #f ])
      (expect [<= 10 (length (cp:get-active-cells pooler)) 40])))
#|                                                                                            ;
(test "InitialInference")
  (let ((pooler (initialize-default-pooler)))
    (cp:compute pooler (range 0 40) #t)      ;; Learn one pattern
    (let ((object-sdr (set (cp:get-active-cells pooler))))
      (cp:compute pooler (range 0 40) #t)    ;; Form internal distal connections
      (cp:reset pooler)
      (cp:compute pooler (range 0 40) #f)    ;; Inferring on same pattern should lead to same result
      (expect [(set (cp:get-active-cells pooler)) object-sdr])
      (cp:compute pooler (list) #f)          ;; Inferring with no inputs should maintain same pattern
      (expect [(set (cp:get-active-cells pooler)) object-sdr])))
                                                                                            ;
(test "ShortInferenceSequence")
  (let ((pooler (initialize-default-pooler)))
    (cp:compute pooler (range 0 40) #t)
    (let ((object1-sdr (set (cp:get-active-cells pooler))))
      (cp:compute pooler (range 100 140) #t)
      (expect [(set (cp:get-active-cells pooler)) object1-sdr])
      (cp:reset pooler)
      (cp:compute pooler (range 1000 1040) #t)
      (let ((object2-sdr (set (cp:get-active-cells pooler))))
        (cp:compute pooler (range 1100 1140) #t)
        (expect [(set (cp:get-active-cells pooler)) object2-sdr])
        (cp:reset pooler)
        (cp:compute pooler (range 100 140) #f)
        (expect [(set (cp:get-active-cells pooler)) object1-sdr])
        (cp:compute pooler (list) #f)
        (expect [(set (cp:get-active-cells pooler)) object1-sdr])
        (cp:reset pooler)
        (cp:compute pooler (range 0 40) #f)
        (expect [(set (cp:get-active-cells pooler)) object1-sdr])
        (cp:reset pooler)
        (cp:compute pooler (range 1100 1140) #f)
        (expect [(set (cp:get-active-cells pooler)) object2-sdr])
        (cp:compute pooler (list) #f)
        (expect [(set (cp:get-active-cells pooler)) object2-sdr])
        (cp:reset pooler)
        (cp:compute pooler (range 1000 1040) #f)
        (expect [(set (cp:get-active-cells pooler)) object2-sdr]))))
                                                                                            ;
(test "ProximalLearning_SampleSize")
  (let ((pooler
          (cp:make-cp `(
            [input-width    . ,(* 2048 8)]
            [initial-proximal-permanence  . ,(perm 0.60)]
            [connected-proximal-permanence . ,(perm 0.50)]
            [sample-size-proximal          . 10]
            [syn-perm-proximal-dec         . 0])))
        (feedforward-input-1 (range 0 10)))
    (cp:compute pooler feedforward-input-1 #t)
    (for-each
      (lambda (cell)
        (expect [(cp:test:number-of-proximal-synapses     pooler (list cell)) 10]
                [(cp:test:num-connected-proximal-synapses pooler (list cell)) 10])
        (let* ( 
            (synapses
              (vector->list (vector-ref (cp:test:cp-proximal-permanences pooler) cell)))
            (presynaptic-cells (map syn-prex synapses))
            (permanences       (map syn-perm synapses)))
          (expect [(set presynaptic-cells) (set feedforward-input-1)])
          (for-each
            (lambda (p)
              (expect [p (perm 0.60)]))
            permanences)))
      (cp:get-active-cells pooler))
    (cp:compute pooler (range 10 20) #t)
    (for-each
      (lambda (cell)
        (expect [(cp:test:number-of-proximal-synapses     pooler (list cell)) 20]))
      (cp:get-active-cells pooler))
    )                                                                                        ;
|#
  ;; flush test failures (replace failures with #t to confirm all tests)
  (when #t ;failures
    (newline)
    (display successes)
    (display " tests passed")
    (display tests)
    (newline))
  (flush-output-port (current-output-port))

