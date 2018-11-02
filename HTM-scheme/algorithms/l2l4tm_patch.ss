#!chezscheme

;; === HTM-scheme L2L4TM Copyright 2018 Roger Turner. ===
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
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

(library (HTM-scheme HTM-scheme algorithms l2l4tm_patch)
                                                                                            ;
(export
  make-patch
  patch-L2s
  patch-L4s
  patch-TMs
  compute
  reset)
                                                                                            ;
(import
  (except (chezscheme) add1 make-list random reset)
          (HTM-scheme HTM-scheme algorithms htm_prelude)
          (HTM-scheme HTM-scheme algorithms htm_concept)
  (prefix (HTM-scheme HTM-scheme algorithms apical_tiebreak_pair_memory) atpm:)
  (prefix (HTM-scheme HTM-scheme algorithms apical_tiebreak_sequence_memory) atsm:)
  (prefix (HTM-scheme HTM-scheme algorithms column_pooler) cp:))
                                                                                            ;
(define-record-type patch (fields        ;; Patch
  L2s                                    ;; CColVecOf CP   [column_pooler]
  L4s                                    ;; CColVecOf ATPM [apical_tiebreak_pair_memory]
  TMs                                    ;; CColVecOf ATSM [apical_tiebreak_sequence_memory]
  enable-feedback)                       ;; Boolean
  (protocol
    (lambda (new)
      (lambda (ncc is nib ef L2-overrides L4-overrides TM-overrides)
          (new
            (build-vector ncc
              (lambda _ (cp:make-cp   (append L2-overrides
                                              (get-default-L2-params is nib)))))
            (build-vector ncc
              (lambda _ (atpm:make-tm (append L4-overrides
                                              (get-default-L4-params is nib)))))
            (build-vector ncc
              (lambda _ (atsm:make-tm (append TM-overrides
                                              (get-default-TM-params is nib)))))
            ef)))))

(define (get-default-L2-params input-size num-input-bits)
  ;; getDefaultL2Params from l2_l4_inference.py
  (let* ( (sample-size-proximal
            (case num-input-bits
              [(20) 10] [(10) 6] [else (int<- (* num-input-bits .6))]))
          (min-threshold-proximal
            (case num-input-bits
              [(20) 5] [(10) 3] [else (int<- (* sample-size-proximal .6))])))
    `(
      [input-width                   . ,(* input-size 16)]
      [cell-count                    . 4096]
      [sdr-size                      . 40]
      [syn-perm-proximal-inc         . ,(perm 0.1)]
      [syn-perm-proximal-dec         . ,(perm 0.001)]
      [initial-proximal-permanence   . ,(perm 0.6)]
      [min-threshold-proximal        . ,min-threshold-proximal]
      [sample-size-proximal          . ,sample-size-proximal]
      [connected-permanence-proximal . ,(perm 0.5)]
      [syn-perm-distal-inc           . ,(perm 0.1)]
      [syn-perm-distal-dec           . ,(perm 0.001)]
      [initial-distal-permanence     . ,(perm 0.41)]
      [activation-threshold-distal   . 13]
      [sample-size-distal            . 20]
      [connected-permanence-distal   . ,(perm 0.5)])))
                                                                                            ;
(define (get-default-L4-params input-size num-input-bits)
  ;; getDefaultL4Params from l2_l4_inference.py
  (let* ( (sample-size (int<- (* 1.5 num-input-bits)))
          (activation-threshold 
            (case num-input-bits
              [(20) 13] [(10) 8] [else (int<- (* num-input-bits .6))]))
          (min-threshold
            (case num-input-bits
              [(20) 13] [(10) 8] [else activation-threshold])))
    `(
      [column-count                       . ,input-size]
      [cells-per-column                   . 16]
      [initial-permanence                 . ,(perm 0.51)]
      [connected-permanence               . ,(perm 0.6)]
      [permanence-increment               . ,(perm 0.1)]
      [permanence-decrement               . ,(perm 0.02)]
      [min-threshold                      . ,min-threshold]
      [basal-predicted-segment-decrement  . ,(perm 0.0)]
      [apical-predicted-segment-decrement . ,(perm 0.0)]
      [activation-threshold               . ,activation-threshold]
      [reduced-basal-threshold            . ,(int<- (* activation-threshold 0.6))]
      [sample-size                        . ,sample-size])))
                                                                                            ;
(define (get-default-TM-params input-size num-input-bits)
  ;; getDefaultTMParams from combined_sequence_experiment.py
  (let* ( (sample-size (int<- (* 1.5 num-input-bits)))
          (activation-threshold 
            (case num-input-bits
              [(20) 18] [(10) 8] [else (int<- (* num-input-bits .6))]))
          (min-threshold
            (case num-input-bits
              [(20) 18] [(10) 8] [else activation-threshold])))
    `(
      [column-count                       . ,input-size]
      [cells-per-column                   . 16]
      [initial-permanence                 . ,(perm 0.41)]
      [connected-permanence               . ,(perm 0.6)]
      [permanence-increment               . ,(perm 0.1)]
      [permanence-decrement               . ,(perm 0.03)]
      [min-threshold                      . ,min-threshold]
      [basal-predicted-segment-decrement  . ,(perm 0.003)]
      [apical-predicted-segment-decrement . ,(perm 0.0)]
      [reduced-basal-threshold            . ,(int<- (* activation-threshold 0.6))]
      [activation-threshold               . ,activation-threshold]
      [sample-size                        . ,sample-size])))
                                                                                            ;
(define (compute patch features locations learn)
  ;; run one time step of a patch with multiple cortical columns; interaction between
  ;; cols is only L2 lateral inputs which are buffered, so can run cols in any order
  (let ((num-columns (vector-length features))
        (columnx (make-thread-parameter 0))
        (mutex   (make-mutex))
        (finished (make-condition))
        (done     0))
    (define (column-thread)                     ;; each cc
      (define (compute-column L2 L4 TM feature location)
        (atpm:compute L4
          feature                        ;; feedforward
          location                       ;; basal input
          (if (patch-enable-feedback patch)
            (cp:get-active-cells L2)     ;; apical input
            '() )
          location                       ;; basal growth candidates
          '()                            ;; apical growth candidates
          learn)
        (cp:compute L2 
          (atpm:get-active-cells L4)     ;; feedforward
          (vector-fold-left              ;; lateral inputs from other cols
            (lambda (lateral-inputs other-L2)
              (if (eq? L2 other-L2)  lateral-inputs
                (cons (cp:get-prev-active-cells other-L2) lateral-inputs)))
            '()
            (patch-L2s patch))
          (atpm:get-predicted-cells L4)  ;; feedforward growth candidates
          learn
          '())                           ;; predicted input
        (atsm:compute TM
          feature                        ;; feedforward
          learn))
      (compute-column
        (vector-ref (patch-L2s patch) (columnx))
        (vector-ref (patch-L4s patch) (columnx))
        (vector-ref (patch-TMs patch) (columnx))
        (vector-ref features  (columnx))
        (vector-ref locations (columnx)))
      (with-mutex mutex
        (set! done (add1 done))
        (when (= done num-columns)
          (condition-signal finished))))
    (do ((threadx 0 (add1 threadx))) ((= threadx num-columns))
      (columnx threadx)
      (fork-thread column-thread))
    (with-mutex mutex
      (condition-wait finished mutex))))
                                                                                            ;
(define (reset patch)
  (vector-for-each
    (lambda (L2 L4 TM)
      (cp:reset   L2)
      (atpm:reset L4)
      (atsm:reset TM))
    (patch-L2s patch) (patch-L4s patch) (patch-TMs patch)))
                                                                                            ;
(define (get-TM-predicted-cols patch)
  (vector-map
    (lambda (TM)
      (map
        (lambda (cellx) (fxdiv cellx 16))
        (atsm:get-predicted-cells TM)))
    (patch-TMs patch)))
                                                                                            ;
(define (get-L2-representations patch)   ;; Patch -> (CColVecOf SDR)
  (vector-map cp:get-active-cells (patch-L2s patch)))
                                                                                            ;
)
