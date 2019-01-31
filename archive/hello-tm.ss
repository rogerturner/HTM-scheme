#!r6rs

;; ========= HTM-scheme Hello-TM example Copyright 2017 Roger Turner. =========
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Based on hello_tm.py which is part of the Numenta Platform for        ;;
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

  ;; Translated from NuPIC hello_tm.py, see comments there for more info.

(import 
  (rnrs)
  (HTM-scheme archive htm_prelude)
  (HTM-scheme archive temporal_memory))

(define (hello-tm)
  (display "See nupic/examples/tm/hello_tm.py") (newline)
  (let* ( (format-row (lambda (x) (vector-fold-left     ;; (vectorof Number) -> String
                                    (lambda (s xc) (string-append s  
                                        (if (zero? (mod (string-length s) 11)) " " "")
                                        (number->string xc)))
                                    "" (vector-take 100 x))))
          (tm (tm:construct '(50) 2                             ;; Step 1: create Temporal Pooler instance with appropriate parameters
                `[initial-permanence    . ,(tm:permanence 0.5)]
                `[connected-permanence  . ,(tm:permanence 0.5)]
                `[min-threshold         . 8]
                `[max-new-synapse-count . 20]
                `[permanence-increment  . ,(tm:permanence 0.1)]
                `[permanence-decrement  . ,(tm:permanence 0.0)]
                `[activation-threshold  . 8]))
          (x (build-vector 5                                    ;; Step 2: create input vectors to feed to the temporal memory. Each input vector
                (lambda (xx)                                    ;; must be numberOfCols wide. Here we create a simple sequence of 5 vectors
                  (build-vector 50                              ;; representing the sequence A -> B -> C -> D -> E
                    (lambda (colx)
                      (if (= (div colx 10) xx) 1 0))))))
          (nzindices (lambda (vec) 
                        (let loop ((i 0) (l '()))
                          (if (= i (vector-length vec)) (reverse l)
                            (loop (add1 i) (if (zero? (vector-ref vec i)) l (cons i l))))))))
                            
    (do ((i 0 (add1 i))) ((= i 10))                             ;; Step 3: send this simple sequence to the temporal memory for learning
      (do ((j 0 (add1 j))) ((= j 5))                            ;; We repeat the sequence 10 times
        (let ((active-columns (nzindices (vector-ref x j))))
          (tm:compute tm active-columns #t)))
      (tm:reset tm))
    (do ((j 0 (add1 j))) ((= j 5))                              ;; Step 4: send the same sequence of vectors and look at predictions made by
      (for-each display `(                                      ;; temporal memory
        "--------" ,(string-ref "ABCDE" j) "--------" #\newline 
        "Raw input vector: " ,(format-row (vector-ref x j)) #\newline))
        (let ((active-columns (nzindices (vector-ref x j))))
        (tm:compute tm active-columns #f))
      (let ((act-col-state (build-vector 50 (lambda (cx)
                               (if (member cx (tm:get-active-cols tm)) 1 0)))))
        (for-each display `("Active columns:   " ,(format-row act-col-state) #\newline)))
      (let ((pred-col-state (build-vector 50 (lambda (cx)
                                (if (member cx (tm:get-predictive-cols tm)) 1 0)))))
        (for-each display `("Predicted columns:" ,(format-row pred-col-state) #\newline))))))
