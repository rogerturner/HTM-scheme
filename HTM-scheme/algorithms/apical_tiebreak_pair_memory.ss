;; === HTM-scheme Apical Tiebreak Pair Memory Copyright 2017 Roger Turner. ===
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

#!chezscheme
                                                                                            ;
(optimize-level 3)
                                                                                            ;
(library (apical_tiebreak_pair_memory)
                                                                                            ;
(export
  tm:permanence
  tm:synapse
  tm:segment
  (rename
    (tm:tm                      atpm:tm)
    (tm:construct               atpm:construct)
    (tm:reset                   atpm:reset)
    (tm:get-active-cells        atpm:get-active-cells)
    (compute                    atpm:compute)
    (get-predicted-cells        atpm:get-predicted-cells)
    (get-basal-predicted-cells  atpm:get-basal-predicted-cells)
    (get-apical-predicted-cells atpm:get-apical-predicted-cells)))
                                                                                            ;
(import 
  (except (chezscheme) add1 make-list random reset)
  (htm_prelude)
  (apical_tiebreak_temporal_memory))
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
    (tm:depolarize-cells tm basal-input apical-input learn)
    (tm:activate-cells tm active-columns basal-input apical-input
      basal-growth-candidates apical-growth-candidates learn)))
                                                                                            ;
(define (get-predicted-cells tm)         ;; TM -> {CellX}
  ;; Cells that were predicted for this timestep
  (tm-predicted-cells tm))
                                                                                            ;
(define (get-basal-predicted-cells tm)   ;; TM -> {CellX}
  ;; Cells with active basal segments
  (unique (map-segments-to-cells (tm-active-basal-segments tm))))
                                                                                            ;
(define (get-apical-predicted-cells tm)  ;; TM -> {CellX}
  ;; Cells with active apical segments
  (unique (map-segments-to-cells (tm-active-apical-segments tm))))
                                                                                            ;
)