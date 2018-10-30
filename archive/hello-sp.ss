;; ========= HTM-scheme Hello-SP example Copyright 2017 Roger Turner. =========
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Based on hello_sp.py which is part of the Numenta Platform for        ;;
  ;; Intelligent Computing (NuPIC) Copyright (C) 2013-2016, Numenta, Inc.  ;;
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

  ;; Translated from NuPIC hello_sp.py, see comments there for more info.

(library-directories "../HTM-scheme/algorithms/")

(import 
  (rnrs)
  (htm_prelude)
  (spatial_pooler))

(define (random-bits size)
  (let ((w (expt 2 32)))
    (do ((i 1 (add1 i)) (r (random w) (+ (* w r) (random w)))) ((= i (div size 32)) r))))
                                                                                            ;
(define (percent-overlap x1 x2)
  (let ((min-x1x2 (min (bitwise-bit-count x1) (bitwise-bit-count x2))))
    (if (positive? min-x1x2)
        (* 100 (/ (bitwise-bit-count (bitwise-and x1 x2)) min-x1x2))
        0)))
                                                                                            ;
(define (percent->string x)
  (string-append (if (< x 10) "  " (if (< x 100) " " "")) 
                 (number->string (inexact (/ (int<- (* 10 x)) 10))) "%"))
                                                                                            ;
(define (corrupt-vector vector size noise-level)
  (let* ( (bit-xs (vector-sample (build-vector size id) (int<- (* noise-level size))))
          (flip   (vector->bitwise bit-xs)))
    (bitwise-xor vector flip)))
                                                                                            ;
(define (hello-sp) 
  (display "See nupic/examples/sp/hello_sp.py") (newline)
  (letrec* ( 
    (inp-dims '(32 32))
    (col-dims '(64 64))
    (ni   (apply * inp-dims))
    (sp   (make-sp* inp-dims col-dims '[global-inhibition . #t]
          `[potential-radius                . ,ni]
          `[num-active-columns-per-inh-area . ,(int<- (* 0.02 (apply * col-dims)))]
          `[syn-perm-active-inc             . ,(perm<- 0.01)]
          `[syn-perm-inactive-dec           . ,(perm<- 0.008)]))
    (create-input (lambda () (random-bits ni)))
    (run          (lambda (description input learn prev-cols)
                    (let* ( (cols  (list-sort < (compute sp input learn)))
                            (score (percent-overlap (list->bitwise cols) (list->bitwise prev-cols))))
                      (for-each display `(,description ": " ,(percent->string score) " ("))
                      (for-each (lambda (x) (display x) (display " ")) (take 15 cols)) 
                      (display "... ") (display (car (reverse cols))) (display ")") (newline)
                      cols)))
    (input1 (create-input))
    (input2 (create-input))
    (input3 (create-input))
    (input4 (corrupt-vector input3 ni 0.1))
    (input5 (corrupt-vector input4 ni 0.2)))
  (let* (
      (cols1 (run "Random 1" input1 #f '()))
      (cols2 (run "Random 2" input2 #f cols1))
      (cols3 (run "Random 3" input3 #t cols2))
      (cols4 (run "Repeat 3" input3 #t cols3)))
    (run "Noise .1" input4 #t cols4)
    (run "Noise .3" input5 #t cols4)))
  'ok)
