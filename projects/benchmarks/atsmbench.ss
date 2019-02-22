;; HTM-scheme atsmbench.ss Copyright 2019 Roger Turner.
 ; This program is released under the GNU Affero General Public License version 3
 ; http://www.gnu.org/licenses.  Derived from marty1885/tiny-htm/examples/tmbench.cpp
  
(import
  (except (chezscheme) add1 make-list reset)
  (except (HTM-scheme HTM-scheme algorithms htm_prelude) random)
  (prefix (HTM-scheme HTM-scheme algorithms apical_tiebreak_sequence_memory) atsm:))

(define (benchmark-temporal-memory       ;; Nat (Vectorof SDR) Nat -> Number Nat Nat
          input-size sdrs num-cc)
  ;; make num-cc sequence memories and run compute on each for each element of sdrs
  ;; produce total thread time and segment and synapse counts
  (let ((time-used   (make-time 'time-monotonic 0 0))
        (thread-colx (make-thread-parameter #f))
        (mutex       (make-mutex))
        (finished    (make-condition))
        (done        0)
        (tms  (build-vector num-cc (lambda _
                (atsm:make-tm `([column-count     . ,input-size]
                                [cells-per-column . 16]))))))
    (vector-for-each (lambda (sdr)
        (define (thread)
          (let ((t0 (current-time 'time-thread)))
            (atsm:compute (vector-ref tms (thread-colx)) sdr #;learn #t)
            (with-mutex mutex
              (set! time-used (add-duration time-used 
                (time-difference (current-time 'time-thread) t0)))
              (set! done (add1 done))
              (when (= done num-cc) (condition-signal finished)))))
        (set! done 0)
        (do ((colx 0 (add1 colx))) ((= colx num-cc))
          (thread-colx colx)
          (fork-thread thread))
        (when (< done num-cc)
          (with-mutex mutex (condition-wait finished mutex))))
      sdrs)
      (values
        (+ (* (time-second time-used) 1000)  ;; total thread time -> milliseconds
              (/ (time-nanosecond time-used) 1000000))
        (vector-fold-left (lambda (acc tm)
            (+ acc (atsm:get-n-segments-created tm)))
          0 tms)
        (vector-fold-left (lambda (acc tm)
            (+ acc (atsm:get-n-synapses-created tm)))
          0 tms))))

(define (encode value min-val max-val    ;; Number Number Number Nat Nat -> SDR
          sdr-length sdr-bits)
  ;; produce SDR encoding of ~length ~bits of value in range min-val - max-val
  (let* ( (encode-space (/ (- sdr-length sdr-bits) (- max-val min-val)))
          (start (int<- (* encode-space (- value min-val)))))
    (build-list sdr-bits (lambda (i) (+ i start)))))

(define (generate-random-data            ;; Nat Nat Number -> (Vectorof SDR)
          input-size num-data sparsity)
  ;; produce a vector with num-data random SDRs of specified input-size and sparsity
  (build-vector num-data (lambda _ 
      (encode (random 1.0) 0 1 input-size 
              (max 10 (int<- (* input-size sparsity)))))))

(define (main)
  (let* ( (args       (append (command-line-arguments) '("" "" "" "")))
          (num-cc     (or (string->number (car args))       4))
          (input-size (or (string->number (cadr args))   1024))
          (sparsity   (or (string->number (caddr args))  0.02))
          (num-data   (or (string->number (cadddr args)) 1000))
          (input-data (generate-random-data input-size num-data sparsity)))
    (for-each display `("HTM-scheme Sequence Memory: " 
                        ,num-cc " columns, " ,input-size " minicols/col, "
                        "sparsity " ,sparsity "\n"))
    (let-values (((t segments synapses)
                  (benchmark-temporal-memory input-size input-data num-cc)))
      (for-each display `(,(inexact (/ (round (* 10 (/ t num-data))) 10)) "ms/compute, "
                          ,segments " segments, " ,synapses " synapses\n")))))
                        
(time (main))
