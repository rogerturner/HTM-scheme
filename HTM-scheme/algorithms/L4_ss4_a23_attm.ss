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
  
Layer4 consists of location and sequence sub-layers (apical tiebreak temporal memories)

               ^    ^         v
               |    |         |
               ^    ^         v       
            +--ac---pc---+-+--ai--+   
      +---->bi   L4loc   |*|      |   
      |     +------------|*|------+               ac = active cells
      |                  |*|                      ai = apical input
      +--------+         |*|                      bi = basal input
      |        ^         |*|                      pc = predicted cells
      |     +--ac--------|*|------+
      +---->bi   L4seq   |*|      |
      |     +------------+-+------+ 
      |                   ^ \
      |                   |  \
      ^                   ^   \
    location        feature   minicolumn spans L4loc and L4seq

    |#

(library (HTM-scheme HTM-scheme algorithms L4_ss4_a23_attm)
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
(define-enumeration L4-populations (loc seq) L4)
                                                                                            ;
(define-record-type l4 (fields           ;; Layer
  loc                                    ;; ATTM ss4/p4, axon -> L2/3
  seq)                                   ;; ATTM ss4, axon -> L4
  (protocol
    (lambda (new)
      (lambda (bis nib ais loc-overrides seq-overrides)
        (new (attm:make-tm (append `([basal-input-size  . ,bis] [apical-input-size . ,ais])
                                    loc-overrides
                                    (get-default-loc-params bis nib)))
             (attm:make-tm (append `([basal-input-size  . ,bis] [apical-input-size . ,ais])
                                    seq-overrides
                                    (get-default-seq-params bis nib))))))))

(define (union l1 l2)
  (if (integer? l1)
    (+ l1 l2)
    (unique! fx=? (union1d l1 l2))))

(define (compute l feature location      ;; Layer SDR SDR SDR Boolean Populations ->
          apical-input learn pop)
  ;; step layer l with sub-layers pop
  (let* ( (loc                     (enum-set-member? 'loc pop))
          (seq                     (enum-set-member? 'seq pop))
          (basal-input
            (union location (attm:get-active-cells (l4-seq l))))
          (basal-growth-candidates (attm:get-winner-cells (l4-seq l))))
    (when loc (attm:depolarize-cells (l4-loc l) basal-input apical-input learn))
    (when seq (attm:depolarize-cells (l4-seq l) basal-input '()          learn))
    (let* (
        (loc-predicted-cells  (if loc (attm:get-predicted-cells (l4-loc l))  '()))
        (loc-predicted-cols   (attm:cols-from-cells (l4-loc l) loc-predicted-cells))
        (seq-predicted-cells  (if seq (attm:get-predicted-cells (l4-seq l))  '()))
        (seq-predicted-cols   (attm:cols-from-cells (l4-seq l) seq-predicted-cells))
        (bursting-columns     (setdiff1d feature (union loc-predicted-cols seq-predicted-cols))))
      (when loc (attm:activate-cells (l4-loc l)
                  feature                ;; feedforward-input
                  basal-input            ;; basal-reinforce-candidates
                  apical-input           ;; apical-reinforce-candidates
                  location #;basal-growth-candidates               ;; basal-growth-candidates
                  apical-input           ;; apical-growth-candidates
                  learn bursting-columns))
      (when seq (attm:activate-cells (l4-seq l)
                  feature                ;; feedforward-input
                  basal-input            ;; basal-reinforce-candidates
                  '()                    ;; apical-reinforce-candidates
                  basal-growth-candidates ;; basal-growth-candidates
                  '()                    ;; apical-growth-candidates
                  learn bursting-columns)))))
                                                                                            ;
(define (reset l)
  (attm:reset (l4-loc l))
  (attm:reset (l4-seq l)))
                                                                                            ;
(define (reset-seq l)
  (attm:reset (l4-seq l)))
                                                                                            ;
(define (get f l pop)
  ;; access f cells of sub-population pop of layer
  (let ((loc  (enum-set-member? 'loc pop))
        (seq  (enum-set-member? 'seq pop)))
    (assert (or loc seq))
    (cond
      [(and loc seq) (union (f (l4-loc l)) (f (l4-seq l))) ]
      [ loc          (f (l4-loc l)) ]
      [ seq          (f (l4-seq l)) ] )))
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