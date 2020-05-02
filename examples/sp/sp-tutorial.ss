;; ====== HTM-scheme SP-Tutorial example Copyright 2017 Roger Turner. ======
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Based on sp_tutorial.py which is part of the Numenta Platform for     ;;
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

  ;; Translated from NuPIC sp_tutorial.py, see comments there for more info.

(import
  (rnrs)
  (HTM-scheme HTM-scheme algorithms htm_prelude)
  (HTM-scheme HTM-scheme algorithms htm_concept)
  (prefix (HTM-scheme HTM-scheme algorithms spatial_pooler) sp:))

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
(define (sp-tutorial)
  (display "See nupic/examples/sp/sp_tutorial.py") (newline)
  (random-seed! 42)
  (let* ( 
      (input-dimensions  '(1024))
      (column-dimensions '(2048))
      (input-size    (apply * input-dimensions))
      (column-number (apply * column-dimensions))
      (create-input  (lambda _ (random-bits input-size)))
      (input-array   (create-input 0))
      (sp (sp:make-sp `(
              [input-dimensions                . ,input-dimensions]
              [column-dimensions               . ,column-dimensions]
              [potential-radius                . ,(int<- (* 0.5 input-size))]
              [num-active-columns-per-inh-area . ,(int<- (* 0.02 column-number))]
              [global-inhibition               . #t]
              [syn-perm-active-inc             . ,(perm 0.01)]
              [syn-perm-inactive-dec           . ,(perm 0.008)])))
      (active-cols   (sp:compute sp input-array #f))
      (all-counts    (sp:calculate-overlap sp input-array))
      (active-counts (vector-refs all-counts (list->vector active-cols)))
      (mean          (lambda (vec) (int<- (vector-average vec))))
      (for-each-noise-level
        (lambda (description proc)
          (display description) (newline)
          (do ((i 0 (add1 i))) ((= i 11) )
            (let ((noise-level (/ i 10.0)))
              (for-each display
                `(,noise-level "  " ,(percent->string (proc noise-level)) #\newline)))))))
    (for-each display 
      `("Figure 1 - overlap count means" #\newline
        "all cols: " ,(mean all-counts) #\newline
        "  active: " ,(mean active-counts) #\newline))
    (for-each-noise-level "Figure 2 - noise:overlap linear" (lambda (nl)
      (percent-overlap input-array (corrupt-vector input-array input-size nl))))
    (for-each-noise-level "Figure 3 - without training" (lambda (nl)
      (percent-overlap (list->bitwise (sp:compute sp input-array #f))
                       (list->bitwise (sp:compute sp 
                         (corrupt-vector input-array input-size nl) #f)))))
    (let* ( (num-examples 10)
            (input-vectors (build-vector num-examples create-input))
            (epochs 30))
      (do ((epoch 0 (add1 epoch))) ((= epoch epochs))
        (do ((i 0 (add1 i))) ((= i num-examples))
          (sp:compute sp (vector-ref input-vectors i) #t)))
      (for-each-noise-level "Figure 4 - with training: sigmoid" (lambda (nl)
        (percent-overlap (list->bitwise (sp:compute sp (vector-ref input-vectors 0) #f))
                         (list->bitwise (sp:compute sp 
                           (corrupt-vector (vector-ref input-vectors 0) input-size nl) #f))))))))
