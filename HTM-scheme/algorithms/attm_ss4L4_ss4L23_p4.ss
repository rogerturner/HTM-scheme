#!chezscheme

;; === HTM-scheme Layer4 algorithm Copyright 2019 Roger Turner. ===
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

  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.
  #|
  
Layer4 consists of ss4L4, ss4L23, and p4 populations (apical tiebreak temporal memories)

               ^      ^       v
               |      |       |
               ^      ^       v
            +--ac-----pc-+-+--ai--+
      +---->bi |  p4  |  |*|      |
      |     +--|------|--|*|------+               ac = active cells
      |        ^      ^  |*|                      (output is union
      |     +--ac-----pc-+*+------+                of p4 and ss4L4)
      +---->bi   ss4L23  |*|      |
      |     +------------|*|------+               ai = apical input
      |                  |*|
      +--------+         |*|                      bi = basal input
      |        ^         |*|                      (union of location
      |     +--ac--------|*|------+                and sss4L4's ac
      +---->bi   ss4L4   |*|      |                on previous step)
      |     +------------+-+------+
      |                   ^ \                     pc = predicted cells
      |                   |  \
      ^                   ^   \
    location        feature   inhibition applies to all cells in minicolumn

    |#

(library (HTM-scheme HTM-scheme algorithms attm_ss4L4_ss4L23_p4)
                                                                                            ;
(export
  L4
  make-l4
  compute
  reset
  reset-seq
  get-active-cells
  get-predicted-cells
  get-predicted-active-cells
  get-n-segments-created
  get-n-synapses-created)
                                                                                            ;
(import
  (except (chezscheme) add1 make-list random reset)
          (HTM-scheme HTM-scheme algorithms htm_prelude)
          (HTM-scheme HTM-scheme algorithms htm_concept)
  (prefix (HTM-scheme HTM-scheme algorithms apical_tiebreak_temporal_memory) attm:))
                                                                                            ;
(define-enumeration L4-populations (ss4L4 ss4L23 p4) L4)
                                                                                            ;
(define-record-type l4 (fields           ;; Layer
  ss4L4                                  ;; ATTM ss4, axon -> L4
  ss4L23                                 ;; ATTM ss4, axon -> L2/3
  p4)                                    ;; ATTM p4
  (protocol
    (lambda (new)
      (lambda (bis nib ais loc-overrides seq-overrides)
        (let ((input-sizes `([basal-input-size  . ,bis] [apical-input-size . ,ais])))
          (new 
            (attm:make-tm (append input-sizes seq-overrides (get-default-seq-params bis nib)))
            (attm:make-tm (append input-sizes loc-overrides (get-default-loc-params bis nib)))
            (attm:make-tm (append input-sizes loc-overrides (get-default-loc-params bis nib)))))))))
                                                                                            ;
(define (compute l feature location      ;; Layer SDR SDR SDR Boolean Populations -> SDR
          apical-input learn pop)
  ;; step layer l with sub-layers pop; return active and predicted L2 projecting cells
  (let* ( (prev-ss4L4 (attm:get-active-cells (l4-ss4L4 l)))
          (basal-input (union1d prev-ss4L4 location)))
    (when (enum-set-member? 'ss4L4 pop) 
      (attm:depolarize-cells (l4-ss4L4 l) basal-input '() learn))
    (when (enum-set-member? 'ss4L23 pop) 
      (attm:depolarize-cells (l4-ss4L23 l) basal-input '() learn))
    (when (enum-set-member? 'p4 pop)
      (attm:depolarize-cells (l4-p4 l) basal-input apical-input learn))
    (let* (
        (ss4L4-predicted-cells  (if (enum-set-member? 'ss4L4 pop) (attm:get-predicted-cells (l4-ss4L4 l))  '()))
        (ss4L4-predicted-cols   (attm:cols-from-cells (l4-ss4L4 l) ss4L4-predicted-cells))
        (ss4L23-predicted-cells (if (enum-set-member? 'ss4L23 pop) (attm:get-predicted-cells (l4-ss4L23 l))  '()))
        (ss4L23-predicted-cols  (attm:cols-from-cells (l4-ss4L23 l) ss4L23-predicted-cells))
        (p4-predicted-cells     (if (enum-set-member? 'p4 pop) (attm:get-predicted-cells (l4-p4 l))  '()))
        (p4-predicted-cols      (attm:cols-from-cells (l4-p4 l) p4-predicted-cells))
        (bursting-columns       (setdiff1d feature
                                  (union1d ss4L4-predicted-cols 
                                    (union1d ss4L23-predicted-cols p4-predicted-cols)))))
      (when (enum-set-member? 'ss4L4 pop)
        (let ((ss4L4-winner-cells (append (attm:get-winner-cells (l4-ss4L4 l)) '())))
          (attm:activate-cells (l4-ss4L4 l)
            feature                              ;; feedforward-input
            basal-input                          ;; basal-reinforce-candidates
            '()                                  ;; apical-reinforce-candidates
            ss4L4-winner-cells   ;; basal-growth-candidates
            '()                                  ;; apical-growth-candidates
            learn 
            bursting-columns)))
      (when (enum-set-member? 'ss4L23 pop)
        (let ((ss4L23-winner-cells (append (attm:get-winner-cells (l4-ss4L23 l)) '())))
          (attm:activate-cells (l4-ss4L23 l)
            feature                              ;; feedforward-input
            basal-input                          ;; basal-reinforce-candidates
            '()                                  ;; apical-reinforce-candidates
            (union1d ss4L23-winner-cells location) ;; basal-growth-candidates
            '()                                  ;; apical-growth-candidates
            learn 
            bursting-columns)))
      (when (enum-set-member? 'p4 pop)
        (let ((p4-winner-cells (append (attm:get-winner-cells (l4-p4 l)) '())))
          (attm:activate-cells (l4-p4 l)
            feature                              ;; feedforward-input
            basal-input                          ;; basal-reinforce-candidates
            apical-input                         ;; apical-reinforce-candidates
            (union1d p4-winner-cells location)   ;; basal-growth-candidates
            apical-input                         ;; apical-growth-candidates
            learn 
            bursting-columns)))
      (values
        (union1d
          (attm:get-active-cells (l4-ss4L23 l))
          (attm:get-active-cells (l4-p4 l)))
        (union1d ss4L23-predicted-cells p4-predicted-cells)))))
                                                                                              ;
(define (reset l)
  (attm:reset (l4-ss4L4 l))
  (attm:reset (l4-ss4L23 l))
  (attm:reset (l4-p4 l)))
                                                                                            ;
(define (reset-seq l)
  (attm:reset (l4-ss4L4 l)))
                                                                                            ;
(define (get f l pop)
  ;; access f cells of sub-population pop of layer
  (cond
    [ (enum-set-member? 'ss4L4 pop)  (f (l4-ss4L4  l)) ]
    [ (enum-set-member? 'ss4L23 pop) (f (l4-ss4L23 l)) ]
    [ (enum-set-member? 'p4 pop)     (f (l4-p4     l)) ] ))
                                                                                            ;
(define (get-active-cells l pop)
  (get attm:get-active-cells l pop))
                                                                                            ;
(define (get-predicted-cells l pop)
  (get attm:get-predicted-cells l pop))
                                                                                            ;
(define (get-predicted-active-cells l pop)
  (get attm:get-predicted-active-cells l pop))
                                                                                            ;
(define (get-n-segments-created l pop)
  (get attm:get-n-segments-created l pop))
                                                                                            ;
(define (get-n-synapses-created l pop)
  (get attm:get-n-synapses-created l pop))
                                                                                            ;
(define (get-default-loc-params input-size num-input-bits)
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
(define (get-default-seq-params input-size num-input-bits)
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
)