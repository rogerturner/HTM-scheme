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
  (mutable L2-prev-actives)              ;; CColVecOf (listof CellX)
  enable-feedback)                       ;; Boolean
  (protocol
    (lambda (new)
      (lambda (ncc is nib ef L2-overrides L4-overrides TM-overrides)
        (let* ( (L2s (build-vector ncc (lambda _ 
                        (cp:make-cp (append L2-overrides (get-default-L2-params is nib))))))
                (L4s (build-vector ncc (lambda _ 
                        (atpm:make-tm (append
                          `([basal-input-size  . ,is]
                            [apical-input-size . ,(cp:number-of-cells (vector-ref L2s 0))])
                          L4-overrides (get-default-L4-params is nib))))))
                (TMs (build-vector ncc (lambda _ 
                        (atsm:make-tm (append
                          `([basal-input-size  . ,is]
                            [apical-input-size . ,(cp:number-of-cells (vector-ref L2s 0))])
                          TM-overrides (get-default-TM-params is nib)))))))
          (new L2s L4s TMs (build-vector ncc (lambda _ '())) ef))))))

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
  ;; run one timestep of a patch with multiple cortical columns
  (let ((num-columns (vector-length (patch-L2s patch)))
        (thread-colx (make-thread-parameter #f))
        (mutex       (make-mutex))
        (finished    (make-condition))
        (done        0))
                                                                                            ;
#;> (define (compute-L4L2 L4 L2 feature location)
      (atpm:compute L4
        feature                        ;; feedforward
        location                       ;; basal input
        (if (patch-enable-feedback patch)
          (cp:get-active-cells L2)     ;; apical input (L2 previous timestep)
          '())
        location                       ;; basal growth candidates
        '()                            ;; apical growth candidates
        learn)
      (cp:compute L2 
        (atpm:get-active-cells L4)     ;; feedforward
        (reverse (vector-fold-left     ;; lateral inputs from other cols
          (lambda (lateral-inputs other-L2 other-prev-active)
            (if (eq? L2 other-L2)  lateral-inputs
              (cons other-prev-active lateral-inputs)))
          '()
          (patch-L2s patch) (patch-L2-prev-actives patch)))
        (atpm:get-predicted-cells L4)  ;; feedforward growth candidates
        learn
        (atpm:get-predicted-cells L4)))                          ;; predicted input
                                                                                            ;
#;> (define (L4L2-thread)
      (compute-L4L2
        (vector-ref (patch-L4s patch) (thread-colx))
        (vector-ref (patch-L2s patch) (thread-colx))
        (vector-ref features  (thread-colx))
        (vector-ref locations (thread-colx)))
      (signal-if-all-done))
                                                                                            ;
#;> (define (TM-thread)
      (atsm:compute
        (vector-ref (patch-TMs patch) (thread-colx))
        (vector-ref features (thread-colx))
        learn)
      (signal-if-all-done))
                                                                                            ;
#;> (define (signal-if-all-done)
      (with-mutex mutex
        (set! done (add1 done))
        (when (= done (* 2 num-columns))
          (condition-signal finished))))
                                                                                            ;
    ;; interaction between cols is L2 lateral inputs which are buffered
    ;; after each timestep, so can run L4+L2/TM for all cols in parallel    
    (do ((colx 0 (add1 colx))) ((= colx num-columns))
      (thread-colx colx)
      (if #f                             ;; #t for single-threaded
        (begin (L4L2-thread) (TM-thread))
        (begin
          (fork-thread L4L2-thread)
          (fork-thread TM-thread))))
    (when (< done (* 2 num-columns))
      (with-mutex mutex
        (condition-wait finished mutex)))
    ;; all columns complete (time t), save activity for t+1 lateral input
    (do ((colx 0 (add1 colx))) ((= colx num-columns))
      (vector-set! (patch-L2-prev-actives patch) colx
        (append                          ;; copy list *elements*
          (cp:get-active-cells (vector-ref (patch-L2s patch) colx))
          '())))))
                                                                                            ;
(define (reset patch)
  (vector-for-each
    (lambda (L2 L4 TM)
      (cp:reset   L2)
      (atpm:reset L4)
      (atsm:reset TM))
    (patch-L2s patch) (patch-L4s patch) (patch-TMs patch))
  (patch-L2-prev-actives-set! patch
    (build-vector (vector-length (patch-L2s patch)) (lambda _ '()))))
                                                                                            ;
)
