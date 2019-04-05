#!chezscheme

;; === HTM-scheme Apical Tiebreak Pair Memory Copyright 2019 Roger Turner. ===
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Based on code from Numenta Platform for Intelligent Computing (NuPIC) ;;
  ;; which is Copyright (C) 2017, Numenta, Inc.                            ;;
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

  ;; Extends the HTM-scheme Apical Tiebreak Temporal Memory library.
  ;; Translated from numenta htmresearch/.../apical_tiebreak_temporal_memory.py --
  ;; see comments there for descriptions of functions and parameters.
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.

(library (HTM-scheme HTM-scheme algorithms apical_tiebreak_pair_memory)
                                                                                            ;
(export
  compute
  get-basal-predicted-cells
  get-apical-predicted-cells
  (rename                            ;; will be renamed to atpm:tm etc on import
    (attm:tm                         tm)
    (attm:make-tm                    make-tm)
    (attm:reset                      reset)
    (attm:number-of-cells            number-of-cells)
    (attm:get-active-cells           get-active-cells)
    (attm:get-predicted-cells        get-predicted-cells)
    (attm:get-predicted-active-cells get-predicted-active-cells)
    (attm:tm-n-segments-created      get-n-segments-created)
    (attm:tm-n-synapses-created      get-n-synapses-created)))
                                                                                            ;
(import 
  (except (chezscheme) add1 make-list random reset)
          (HTM-scheme HTM-scheme algorithms htm_prelude)
          (HTM-scheme HTM-scheme algorithms htm_concept)
  (prefix (HTM-scheme HTM-scheme algorithms apical_tiebreak_temporal_memory) attm:))
                                                                                            ;
(define (compute tm active-columns       ;; TM {ColX} {CellX} {CellX} {CellX} {CellX} Boolean ->
          basal-input apical-input basal-growth-candidates apical-growth-candidates learn)
  ;; Perform one timestep. Use the basal and apical input to form a set of
  ;; predictions, then activate the specified columns, then learn.
  (let ((basal-growth-candidates
          (if (null? basal-growth-candidates)
              basal-input
              basal-growth-candidates))
        (apical-growth-candidates
          (if (null? apical-growth-candidates)
              apical-input
              apical-growth-candidates)))
    (attm:depolarize-cells tm basal-input apical-input learn)
    (attm:activate-cells tm active-columns basal-input apical-input
      basal-growth-candidates apical-growth-candidates learn)))
                                                                                            ;
(define (get-basal-predicted-cells tm)   ;; TM -> {CellX}
  ;; Cells with active basal segments
  (unique fx=? (attm:map-segments-to-cells (attm:get-active-basal-segments tm))))
                                                                                            ;
(define (get-apical-predicted-cells tm)  ;; TM -> {CellX}
  ;; Cells with active apical segments
  (unique fx=? (attm:map-segments-to-cells (attm:get-active-apical-segments tm))))
                                                                                            ;
)
