#!r6rs

;; === HTM-scheme Apical Tiebreak Sequence Memory Copyright 2017 Roger Turner. ===
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

(library (HTM-scheme HTM-scheme algorithms apical_tiebreak_sequence_memory)
                                                                                            ;
(export
  tm
  make-tm
  reset
  compute
  get-predicted-cells
  get-next-basal-predicted-cells
  get-next-apical-predicted-cells
  (rename
    (attm:get-active-cells           get-active-cells)
    (attm:get-predicted-cells        get-next-predicted-cells)
    (attm:get-predicted-active-cells get-predicted-active-cells)))
                                                                                            ;
(import 
  (rnrs)
          (HTM-scheme HTM-scheme algorithms htm_prelude)
          (HTM-scheme HTM-scheme algorithms htm_concept)
  (prefix (HTM-scheme HTM-scheme algorithms apical_tiebreak_temporal_memory) attm:))
                                                                                            ;
(define-record-type tm                   ;; TM
  (parent attm:tm)
  (fields
    (mutable prev-apical-input)             ;; (listof CellX)
    (mutable prev-apical-growth-candidates) ;; (listof CellX)
    (mutable prev-predicted-cells))         ;; (listof CellX)
  (protocol
    (lambda (pargs->new)                 ;; (listof KWarg) -> TM
      (lambda (kwargs)
        (let ((tm (apply (pargs->new kwargs)
                         (key-word-args kwargs atsm-defaults))))
          (attm:tm-basal-input-size-set! tm (* (attm:tm-column-count tm) (attm:tm-cells-per-column tm)))
          tm)))))
                                                                                           ;
(define atsm-defaults                    ;; (listof KWarg)
  `(
    [prev-apical-input                  . ()]
    [prev-apical-growth-candidates      . ()]
    [prev-predicted-cells               . ()]))
    
                                                                                            ;
(define (reset tm)                       ;; TM ->
  (attm:reset tm)
  (tm-prev-apical-input-set! tm             '())
  (tm-prev-apical-growth-candidates-set! tm '())
  (tm-prev-predicted-cells-set! tm          '()))
                                                                                            ;
(define compute                          ;; TM {ColX} {CellX} [{CellX}] Boolean ->
  (case-lambda
  [(tm active-columns apical-input apical-growth-candidates learn)
  ;; Perform one timestep. Activate the specified columns, using the predictions
  ;; from the previous timestep, then learn. Then form a new set of predictions
  ;; using the new active cells and the apicalInput.
    (let ((apical-growth-candidates
            (if (null? apical-growth-candidates)
                apical-input
                apical-growth-candidates)))
      (tm-prev-predicted-cells-set! tm (attm:get-predicted-cells tm))
      (attm:activate-cells tm active-columns (attm:get-active-cells tm) (tm-prev-apical-input tm)
                      (attm:get-winner-cells tm) (tm-prev-apical-growth-candidates tm) learn)
      (attm:depolarize-cells tm (attm:get-active-cells tm) apical-input learn)
      (tm-prev-apical-input-set! tm apical-input)
      (tm-prev-apical-growth-candidates-set! tm apical-growth-candidates))]
  [(tm active-columns apical-input learn)
    (compute tm active-columns apical-input apical-input learn)]
  [(tm active-columns learn)
    (compute tm active-columns '() learn)]))
                                                                                           ;
(define (get-predicted-cells tm)         ;; TM -> {CellX}
  ;; The prediction from the previous timestep
  (tm-prev-predicted-cells tm))
                                                                                            ;
(define (get-next-basal-predicted-cells tm)  ;; TM -> {CellX}
  ;; Cells with active basal segments
  (unique fx=? (attm:map-segments-to-cells (attm:get-active-basal-segments tm))))
                                                                                            ;
(define (get-next-apical-predicted-cells tm) ;; TM -> {CellX}
  ;; Cells with active apical segments
  (unique fx=? (attm:map-segments-to-cells (attm:get-active-apical-segments tm))))

)