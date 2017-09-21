;; ====== HTM-scheme TM-High-Order example Copyright 2017 Roger Turner. ======
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Based on tm_high_order.py which is part of the Numenta Platform for   ;;
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

  ;; Translated from NuPIC tm_high_order.py, see comments there for more info.

(library-directories "../src/")

(import (rnrs)
        (libraries htm-prelude)
        (libraries htm-tm))

(define (cols->string cols)              ;; (listof ColX) -> String [of hex chars]
  (list->string (reverse (string->list (number->string (list->bitwise cols) 16)))))
                                                                                            ;
(define (accuracy current predicted)     ;; (listof ColX) (listof ColX) -> Number
  (let ((num-predicted-cols (length predicted))
        (num-common-cols (bitwise-bit-count (bitwise-and
                                             (list->bitwise current)
                                             (list->bitwise predicted)))))
    (inexact (if (positive? num-predicted-cols)
                (/ num-common-cols num-predicted-cols)
                0))))
                                                                                            ;
(define (corrupt-vector v1 noise-level num-active-cols) ;; InputVec Number ColX -> InputVec
  (let* ( (bit-xs (vector-sample (build-vector 2048 id) 
                                 (int<- (* noise-level num-active-cols))))
          (flip   (vector->bitwise bit-xs)))
    (bitwise-xor v1 flip)))
                                                                                            ;
(define (show-predictions tm sequence)
  (do ((k 0 (add1 k))) ((= k 6))
    (tm:reset tm)
    (tm:compute tm (bitwise->list (vector-ref sequence k)) #f)
    (for-each display `(
                        "--- " ,(string-ref "ABCDXY" k) " ---" #\newline))
    (for-each display `(
                        "   Active cols: " ,(cols->string (tm:get-active-cols tm)) #\newline))
    (for-each display `(
                        "Predicted cols: " ,(cols->string (tm:get-predictive-cols tm)) #\newline))))
                                                                                            ;
(define (train tm sequence time-steps noise-level)
  (let ((accuracies '())
        (prev-average 0))
    (display "Average accuracy: ")
    (do ((t 0 (add1 t))) ((= t time-steps))
      (tm:reset tm)
      (let ((predicted-cols '()))
        (do ((k 0 (add1 k))) ((= k 4))
          (let ((v (corrupt-vector (vector-ref sequence k) noise-level (exact (truncate (* 2048 0.02))))))
            (tm:compute tm (bitwise->list v) #t))
          (if (positive? k)
            (set! accuracies (cons (accuracy (tm:get-active-cols tm) predicted-cols) accuracies)))
          (set! predicted-cols (tm:get-predictive-cols tm)))
        (let ((average (/ (round (* 100 (list-average accuracies))) 100)))
          (if (= average prev-average)
              (display ".")
              (begin
                (display average)
                (display " ")
                (set! prev-average average))))))
    (newline)))
                                                                                            ;
(define (tm-high-order)
  (display "See nupic/examples/tm/tm-high-order.py") (newline)
  (let* ( (tm (tm:constructor '(2048) 8
                `[initial-permanence          . ,(tm:permanence 0.21)]
                `[connected-permanence        . ,(tm:permanence 0.5)]
                `[min-threshold               . 10]
                `[max-new-synapse-count       . 20]
                `[permanence-increment        . ,(tm:permanence 0.1)]
                `[permanence-decrement        . ,(tm:permanence 0.1)]
                `[activation-threshold        . 13]
                `[predicted-segment-decrement . ,(tm:permanence 0.004)]))
          (sparsity 0.02)
          (sparse-cols (exact (truncate (* 2048 sparsity))))
          (bits (lambda (b) 
                  (bitwise-copy-bit-field 0 (* b sparse-cols) (* (add1 b) sparse-cols) -1)))
          (sdr-A (bits 0))
          (sdr-B (bits 1))
          (sdr-C (bits 2))
          (sdr-D (bits 3))
          (sdr-X (bits 4))
          (sdr-Y (bits 5))
          (seq1 (vector sdr-A sdr-B sdr-C sdr-D))
          (seq2 (vector sdr-X sdr-B sdr-C sdr-Y))
          (seqT (vector sdr-A sdr-B sdr-C sdr-D sdr-X sdr-Y)))
    
    (display "\nPart 1: A->B->C->D\n")
    (train tm seq1 10 0.0)
    (show-predictions tm seqT)

    (display "\nPart 2: X->B->C->Y\n")
    (train tm seq2 10 0.0)
    (show-predictions tm seqT)
    
    (display "\nPart 3: X->B->C->Y with 30% noise\n")
    (train tm seq2 50 0.3)
    (show-predictions tm seqT)

    (display "\nPart 3 fig 4: X->B->C->Y with 50% noise\n")
    (train tm seq2 50 0.5)
    (show-predictions tm seqT)
    
    (display "\nPart 3 fig 5: X->B->C->Y, no noise\n")
    (train tm seq2 10 0.0)
    (show-predictions tm seqT)
    
    (display "\nPart 3 fig 6: X->B->C->Y with 90% noise\n")
    (train tm seq2 50 0.9)
    (show-predictions tm seqT)
    
    (display "\nPart 3 fig 7: X->B->C->Y, no noise\n")
    (train tm seq2 25 0.0)
    (show-predictions tm seqT)
  
    (display "\nPart 4: mixed ABCD XBCY\n")
    (let ((tm (tm:constructor '(2048) 8
                `[initial-permanence          . ,(tm:permanence 0.21)]
                `[connected-permanence        . ,(tm:permanence 0.5)]
                `[min-threshold               . 10]
                `[max-new-synapse-count       . 20]
                `[permanence-increment        . ,(tm:permanence 0.1)]
                `[permanence-decrement        . ,(tm:permanence 0.1)]
                `[activation-threshold        . 13]
                `[predicted-segment-decrement . ,(tm:permanence 0.004)])))
      (do ((t 0 (add1 t))) ((= t 75))
        (let ((rnd (random 2)))
          (do ((k 0 (add1 k))) ((= k 4))
            (if (zero? rnd)
              (tm:compute tm (bitwise->list (vector-ref seq1 k)) #t)
              (tm:compute tm (bitwise->list (vector-ref seq2 k)) #t)))))
      (show-predictions tm seqT))))
    
