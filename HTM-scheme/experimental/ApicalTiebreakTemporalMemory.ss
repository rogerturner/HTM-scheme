;; === HTM-scheme Apical Tiebreak Temporal Memory Copyright 2017 Roger Turner. ===
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Based on code from Numenta Platform for Intelligent Computing (NuPIC) ;;
  ;; which is Copyright (C) 2013-2016, Numenta, Inc.                       ;;
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

  ;; Extends the HTM-scheme Temporal Memory library htm-tm.ss -- extends
  ;; TM record and uses htm-tm procedures where applicable.
  ;; Translated from NuPIC ApicalTiebreakTemporalMemory.hpp and .cpp, and
  ;; ApicalTiebreakTemporalMemoryTest.cpp -- see comments there for 
  ;; descriptions of functions and parameters.
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.

#!chezscheme

;(optimize-level 3)
                                                                                            ;
(library (libraries htm-attm)
                                                                                            ;
(export
  tm:perm
  tm:prex
  tm:permanence
  tm:synapse
  tm:segment
  tm:get-active-cols
  (rename
    (make-tm              tm:constructor)
    (compute              tm:compute)
    (reset                tm:reset)
    (get-predictive-cols  tm:get-predictive-cols)
    (get-predicted-cells  tm:get-predicted-cells)))
                                                                                            ;
(import 
  (except (chezscheme) add1 make-list random reset)
  (libraries htm-prelude)
  (rename (libraries htm-tm)
    (tm-active-segments             tm-active-basal-segments)
    (tm-active-segments-set!        tm-active-basal-segments-set!)
    (tm-matching-segments           tm-matching-basal-segments)
    (tm-matching-segments-set!      tm-matching-basal-segments-set!)
    (tm-cells                       tm-basal-connections)
    (tm-pre-index                   tm-basal-pre-index)
    (tm-predicted-segment-decrement tm-basal-predicted-segment-decrement)
    ))

;; === Temporal Memory Types ===
                                                                                            ;
;; Boolean, Number, Integer, Fixnum, (listof X), (vectorof X) ... = Scheme types
;; X, Y, Z      = type parameters (arbitrary types in function type specification etc)
;; X Y -> Z     = function with argument types X, Y and result type Z
;; {X}          = abbreviation for (listof X)
;; Nat          = natural number (including zero) (Scheme Fixnum or exact Integer)
;; Fix10k       = Fixnum interpreted as number with 4 decimal places, ie 10000 = 1.0
;; Permanence   = Fix10k [0 - MAXPERM], interpreted as 0.0 - 0.9999
;; Synapse      = Fixnum, interpreted as CellX * 10000 + Permanence 
;; Synapses     = (vectorof Synapse) [optionally fxvector in Chez Scheme]
;; Segment      = Record: CellX, FlatX, and Synapses
;; FlatX        = Nat, segment sequence number
;; CellX        = Nat [0 - (div (expt 2 (fixnum-width)) 10000)] cell index
;; CellVecOf    = Vector indexed by CellX
;; ColX         = Nat [0 - MAXCOL], column index of cell
;; Connections  = (CellVecOf (listof Segment)), basal and apical segments for cells
;; TM           = Record: tm parameters, cells

;; === Parameters, Data, Convenience Functions ===
                                                                                            ;
(define-record-type tm                   ;; TM
  (parent tm:tm)
  (fields
    basal-input-size                     ;; Nat
    apical-input-size                    ;; Nat
    sample-size                          ;; Nat
    ;;basal-predicted-segment-decrement  ;; Permanence [renamed tm predicted-segment-decrement]
    apical-predicted-segment-decrement   ;; Permanence
    (mutable predicted-cells)            ;; (listof CellX)
    (mutable predicted-active-cells)     ;; (listof CellX) 
    ;;(mutable active-basal-segments)    ;; (listof Segment) [renamed tm active-segments]
    ;;(mutable matching-basal-segments)  ;; (listof Segment) [renamed tm matching-segments]
    ;;basal-connections                  ;; (CellVecOf (listof Segment)) [renamed tm cells]
    (mutable basal-potential-overlaps)   ;; (vectorof Nat) [indexed by flatx]
    (mutable active-apical-segments)     ;; (listof Segment)
    (mutable matching-apical-segments)   ;; (listof Segment)
    apical-connections                   ;; Connections
    apical-pre-index                     ;; (CellVecOf {FlatX})
    (mutable apical-potential-overlaps)  ;; (vectorof Nat) [indexed by flatx]
    learn-on-one-cell                    ;; Boolean
    (mutable chosen-cell-for-column)     ;; #f|(ColVecOf CellX|#f)
    (mutable prev-predicted-cells)       ;; (listof CellX)
    (mutable prev-apical-input)          ;; (listof CellX)
    (mutable prev-apical-growth-candidates) ;; (listof CellX)
    )
  (protocol
    (lambda (pargs->new)
      (lambda (cd cpc . args)            ;; (listof Nat) Nat (listof KWarg) -> TM
        (let* ( (num-columns (apply * cd))
                (num-cells   (* num-columns cpc)))
          (apply (apply pargs->new (append (list cd cpc) args)) 
                 (key-word-args
                    (append `(
                      [apical-connections     . ,(make-vector num-cells '())]
                      [apical-pre-index       . ,(make-vector num-cells '())]
                      [chosen-cell-for-column . ,(if (assq 'learn-on-one-cell args)
                                                     (make-vector num-columns #f)
                                                     #f) ])
                      args)
                    attm-defaults)))))))
                                                                                            ;
(define attm-defaults                    ;; (listof KWarg)
  `(
    [basal-input-size                   . 0]
    [apical-input-size                  . 0]
    [sample-size                        . 100]
    [apical-predicted-segment-decrement . ,(tm:permanence 0.004)]
    [predicted-cells                    . ()]
    [predicted-active-cells             . ()]
    [basal-potential-overlaps           . #()]
    [active-apical-segments             . ()]
    [matching-apical-segments           . ()]
    [apical-connections                 . #()]
    [apical-pre-index                   . #()]
    [apical-potential-overlaps          . #()]
    [learn-on-one-cell                  . #f]
    [chosen-cell-for-column             . #f]
    [prev-predicted-cells               . ()]
    [prev-apical-input                  . ()]
    [prev-apical-growth-candidates      . ()]
    ))
                                                                                            ;
(define (cellxs->colx tm cellxs)         ;; TM {CellX} -> ColX
  ;; produce first column index from list, or terminating value if list null
  (if (null? cellxs) (tm-num-columns tm) 
      (cellx->colx tm (car cellxs))))
                                                                                            ;
(define (advance-to-cell tm cellx segs)  ;; TM CellX {Segment} -> {Segment}
  ;; produce next segment list element for given cell
  (let loop ((segs segs))
    (if (null? segs) segs
      (let ((nextseg-cellx (segs->cellx tm (cdr segs))))
        (cond [(fx=? nextseg-cellx cellx) (cdr segs)]
              [(fx<? nextseg-cellx cellx) (loop (cdr segs))]
              [ else segs])))))
                                                                                            ;
(define (advance-to-col tm colx segs)    ;; TM ColX {Segment} -> {Segment}
  ;; produce next segment list element for given column
  (let loop ((segs segs))
    (if (null? segs) segs
      (let ((nextseg-colx (cellx->colx tm (segs->cellx tm (cdr segs)))))
        (cond [(fx=? nextseg-colx colx) (cdr segs)]
              [(fx<? nextseg-colx colx) (loop (cdr segs))]
              [ else segs])))))
                                                                                            ;
(define (advance-cell-to-col             ;; TM ColX {CellX} -> {CellX}
          tm colx cellxs)
  ;; produce next list element for given column
  (let loop ((cellxs cellxs))
    (if (null? cellxs) cellxs
      (let ((next-colx (cellx->colx tm (car cellxs))))
        (cond [(fx=? next-colx colx) cellxs]
              [(fx<? next-colx colx) (loop (cdr cellxs))]
              [ else cellxs])))))

;; === Apical Tiebreak Temporal Memory Algorithm ===
                                                                                            ;
(define (learn-on-cell                   ;; TM Connections CellX {Seg} {Seg} ... ->
          tm connections pre-index cellx active-segments matching-segments 
          active-input potential-overlaps growth-candidates)
  ;; adjust permanences and add new synapses to segments
  (let ((adjust-segment 
          (lambda (segment)
            (adapt-segment tm segment (list->vector active-input) connections)
            (let ((n-grow-desired (fx- (tm-sample-size tm)
                                       (vector-ref potential-overlaps (seg-flatx segment)))))
              (when (fxpositive? n-grow-desired)
                (grow-synapses tm segment n-grow-desired growth-candidates pre-index))))))
    (if (fx=? cellx (segs->cellx tm active-segments))
        ;; Learn on every active segment.
        (let loop ((segs active-segments))
          (when (fx=? cellx (segs->cellx tm segs))
            (adjust-segment (car segs))
            (loop (cdr segs))))
        ;; No active segments. Learn on the best matching segment.
        (if (fx=? cellx (segs->cellx tm matching-segments))
            (let loop ((segs matching-segments) (best-seg '()) (max-overlap (least-fixnum)))
              (if (fx=? cellx (segs->cellx tm segs))
                  (let ((overlap (vector-ref potential-overlaps (seg-flatx (car segs)))))
                    (if (fx>? overlap max-overlap)
                        (loop (cdr segs) segs overlap)
                        (loop (cdr segs) best-seg max-overlap)))
                  (adjust-segment (car best-seg))))
            ;; No matching segments. Grow a new segment and learn on it.
            ;; Don't grow a segment that will never match.
            (let ((n-grow-exact (fxmin (tm-sample-size tm) (length growth-candidates))))
              (when (fxpositive? n-grow-exact)
                (let ((segment (create-segment tm cellx connections)))
                  (grow-synapses tm segment n-grow-exact growth-candidates pre-index))))))))
                                                                                            ;
(define (activate-predicted-column       ;; TM {CellX} {CellX} {CellX} ColX ... -> {CellX} {CellX} {CellX}
          tm active-cells winner-cells predicted-active-cells colx predicted-cells
          active-basal-segs matching-basal-segs active-apical-segs matching-apical-segs
          basal-input basal-growth-candidates
          apical-input apical-growth-candidates learn)
  ;; add predicted cells for column to active/winner/predicted-active lists
  (let loop ( (predicted-cells        predicted-cells)
              (active-basal-segs      active-basal-segs)
              (matching-basal-segs    matching-basal-segs)
              (active-apical-segs     active-apical-segs)
              (matching-apical-segs   matching-apical-segs)
              (active-cells           active-cells) 
              (winner-cells           winner-cells)
              (predicted-active-cells predicted-active-cells))
    (if (or (null? predicted-cells) (not (fx=? (cellx->colx tm (car predicted-cells)) colx)))
      ;; end of predicted cells for this column
      (values active-cells winner-cells predicted-active-cells)
      (let ((cellx (car predicted-cells)))
        (when (fx=? (cellx->colx tm cellx) colx)
          (when learn
            (learn-on-cell tm (tm-basal-connections tm) (tm-basal-pre-index tm) cellx
                           active-basal-segs matching-basal-segs basal-input
                           (tm-basal-potential-overlaps tm) basal-growth-candidates)
            (learn-on-cell tm (tm-apical-connections tm) (tm-apical-pre-index tm) cellx
                           active-apical-segs matching-apical-segs apical-input
                           (tm-apical-potential-overlaps tm) apical-growth-candidates))
          (let* ( (predicted-cells (cdr predicted-cells))
                  (next-cellx (if (null? predicted-cells) (tm-num-cells tm)
                                                          (car predicted-cells))))
            (let ((active-basal-segs    (advance-to-cell tm next-cellx active-basal-segs))
                  (matching-basal-segs  (advance-to-cell tm next-cellx matching-basal-segs))
                  (active-apical-segs   (advance-to-cell tm next-cellx active-apical-segs))
                  (matching-apical-segs (advance-to-cell tm next-cellx matching-apical-segs)))
              (loop predicted-cells active-basal-segs matching-basal-segs active-apical-segs matching-apical-segs
                    (cons cellx active-cells) (cons cellx winner-cells) (cons cellx predicted-active-cells)))))
                    ))))
                                                                                            ;
(define (segment-with-max-overlaps       ;; TM CellX {Segment} (vectorof Nat) -> {Segment}
          tm colx segments overlaps)
  (let loop ((segments segments) (best-seg '()) (max-overlap (least-fixnum)))
    (if (fx=? colx (segs->colx tm segments))
        (let ((overlap (vector-ref overlaps (seg-flatx (car segments)))))
          (if (fx>? overlap max-overlap)
              (loop (cdr segments) segments overlap)
              (loop (cdr segments) best-seg max-overlap)))
        best-seg)))
                                                                                            ;
(define (burst-column                    ;; TM {CellX} {CellX} ColX ... -> {CellX} {CellX}
          tm active-cells winner-cells colx 
          active-basal-segs matching-basal-segs active-apical-segs matching-apical-segs
          basal-input basal-potential-overlaps basal-growth-candidates
          apical-input apical-potential-overlaps apical-growth-candidates learn)
  ;; choose winner cell and learn on it; add all cells in column to active cells
  (let* ( (last             (fx- (fx* (add1 colx) (tm-cells-per-column tm)) 1))
          (cells-for-column (build-list (tm-cells-per-column tm)
                                        (lambda (x) (fx- last x))))
          (chosen-cell      (if (tm-learn-on-one-cell tm)
                                (vector-ref (tm-chosen-cell-for-column tm) colx)
                                #f)))
    (let-values ([(basal-candidates winner-cell)
      (if chosen-cell
          (values matching-basal-segs chosen-cell)
          (if (fx=? colx (segs->colx tm matching-basal-segs))
              (let ((best-basal-segment
                  (segment-with-max-overlaps tm colx matching-basal-segs basal-potential-overlaps)))
                (values best-basal-segment (segs->cellx tm best-basal-segment)))
              (values matching-basal-segs (least-used-cell tm cells-for-column))))])
      (when (tm-learn-on-one-cell tm)
        (vector-set! (tm-chosen-cell-for-column tm) colx winner-cell))
      (when learn
        (learn-on-cell tm (tm-basal-connections tm) (tm-basal-pre-index tm) winner-cell 
                       (advance-to-cell tm winner-cell active-basal-segs)
                       (advance-to-cell tm winner-cell basal-candidates)
                       basal-input basal-potential-overlaps basal-growth-candidates) 
        (learn-on-cell tm (tm-apical-connections tm) (tm-apical-pre-index tm) winner-cell 
                       (advance-to-cell tm winner-cell active-apical-segs)
                       (advance-to-cell tm winner-cell matching-apical-segs)
                       apical-input apical-potential-overlaps apical-growth-candidates))
      (values
        (append cells-for-column active-cells)
        (cons winner-cell winner-cells)))))
                                                                                            ;
(define (activate-cells                  ;; TM {ColX} {CellX} {CellX} {CellX} {CellX} Boolean ->
          tm active-columns basal-reinforce-candidates apical-reinforce-candidates
          basal-growth-candidates apical-growth-candidates learn)
  ;; step through columns in merge of active columns, predicted cells, and matching 
  ;; segments, building active/winner/predicted-active cell lists; this blends in
  ;; the required aspects of iterGroupBy logic from GroupBy.hpp
  ;; active/winner/predicted-active cell lists are current time step,
  ;; tm active/winner/predicted-active are previous time step lists until end
  (let loop ( (active-cols            active-columns)
              (predicted-cells        (tm-predicted-cells tm))
              (active-basal-segs      (tm-active-basal-segments tm)) 
              (matching-basal-segs    (tm-matching-basal-segments tm)) 
              (active-apical-segs     (tm-active-apical-segments tm)) 
              (matching-apical-segs   (tm-matching-apical-segments tm)) 
              (active-cells           '()) 
              (winner-cells           '())
              (predicted-active-cells '()))
    (let ((actcols-colx   (colxs->colx  tm active-cols))
          (predcells-colx (cellxs->colx tm predicted-cells))
          (matbsegs-colx  (segs->colx   tm matching-basal-segs))
          (matasegs-colx  (segs->colx   tm matching-apical-segs)))
      (let ((next (lambda (not-active? actcells wincells predactcells)
                    (let* ( (active-cols 
                              (if (or not-active? (null? active-cols)) active-cols 
                                  (cdr active-cols)))
                            (next-colx (colxs->colx tm active-cols))
                            (matching-basal-segs
                              (if not-active?
                                  (skip-col tm matbsegs-colx matching-basal-segs)
                                  (advance-to-col tm next-colx matching-basal-segs)))
                            (matching-apical-segs
                              (if not-active?
                                  (skip-col tm matasegs-colx matching-apical-segs)
                                  (advance-to-col tm next-colx matching-apical-segs)))
                            (next-colx (fxmin next-colx
                                              (segs->colx  tm matching-basal-segs)
                                              (segs->colx  tm matching-apical-segs))))
                      (loop active-cols
                            (advance-cell-to-col tm next-colx predicted-cells) 
                            (advance-to-col tm next-colx active-basal-segs) 
                            matching-basal-segs 
                            (advance-to-col tm next-colx active-apical-segs) 
                            matching-apical-segs 
                            actcells wincells predactcells))))
            (colx (fxmin actcols-colx predcells-colx matbsegs-colx matasegs-colx)))
        (if (fx=? colx (tm-num-columns tm))
          ;; all active cols and segs handled: set prev active/winner cells for next iteration
          (begin
            (tm-active-cells-set! tm (list->vector (reverse active-cells)))
            (tm-winner-cells-set! tm (reverse winner-cells))
            (tm-predicted-active-cells-set! tm predicted-active-cells))
          (if (fx=? colx actcols-colx)
              (if (fx=? colx predcells-colx)
                  ;; active and predicted
                  (let-values 
                    ([(active-cells winner-cells predicted-active-cells)
                      (activate-predicted-column tm 
                        active-cells winner-cells predicted-active-cells colx predicted-cells
                        active-basal-segs matching-basal-segs active-apical-segs matching-apical-segs
                        basal-reinforce-candidates basal-growth-candidates apical-reinforce-candidates apical-growth-candidates learn)])
                    (next #f active-cells winner-cells predicted-active-cells))
                  ;; else: active but not predicted
                  (let-values 
                    ([(active-cells winner-cells)
                      (burst-column tm active-cells winner-cells colx
                        active-basal-segs matching-basal-segs active-apical-segs matching-apical-segs
                        basal-reinforce-candidates (tm-basal-potential-overlaps tm) basal-growth-candidates 
                        apical-reinforce-candidates (tm-apical-potential-overlaps tm) apical-growth-candidates learn)])
                    (next #f active-cells winner-cells predicted-active-cells)))
            ;; else: not active but colx in predicted|*-segs
            (begin
              (when learn
                (punish-predicted-column tm 
                  (segs->colx tm matching-basal-segs) matching-basal-segs
                  (list->vector basal-reinforce-candidates) (tm-basal-predicted-segment-decrement tm))
                (punish-predicted-column tm 
                  (segs->colx tm matching-apical-segs) matching-apical-segs
                  (list->vector apical-reinforce-candidates) (tm-apical-predicted-segment-decrement tm)))
              (next #t active-cells winner-cells predicted-active-cells))))))))
                                                                                            ;
(define (calculate-predicted-cells       ;; TM {Segment} {Segment} -> {CellX}
          tm basal-segs apical-segs)
  ;; for listed segs, cells with max depolarization score within each column are predicted
  (let each-col ((basal-segs basal-segs) (apical-segs apical-segs) 
                 (colx (fxmin (segs->colx tm basal-segs) (segs->colx tm apical-segs)))
                 (predicted-cells '()))
    (let ((max-depolarization
      (let loop ((basal-segs basal-segs) (apical-segs apical-segs) (max-this-col 0))
        (let* ( (basal-cellx  (segs->cellx tm basal-segs))
                (apical-cellx (segs->cellx tm apical-segs))
                (cellx (fxmin basal-cellx apical-cellx)))
          (if (and (fx<? cellx (tm-num-cells tm))
                   (fx=? (cellx->colx tm cellx) colx))
            (loop 
              (if (fx=? cellx basal-cellx)  (cdr basal-segs)  basal-segs)
              (if (fx=? cellx apical-cellx) (cdr apical-segs) apical-segs)
              (fxmax max-this-col (fx+ (if (fx=? cellx basal-cellx)  2 0)
                                       (if (fx=? cellx apical-cellx) 1 0))))
            max-this-col)))))
      (let loop ((basal-segs basal-segs) (apical-segs apical-segs)
                 (predicted-cells-for-col '()))
        (let* ( (basal-cellx  (segs->cellx tm basal-segs))
                (apical-cellx (segs->cellx tm apical-segs))
                (cellx (fxmin basal-cellx apical-cellx))
                (this-colx (cellx->colx tm cellx)))
          (if (and (fx<? cellx (tm-num-cells tm))
                   (fx=? this-colx colx))
            (loop 
              (if (fx=? cellx basal-cellx)  (cdr basal-segs)  basal-segs)
              (if (fx=? cellx apical-cellx) (cdr apical-segs) apical-segs)
              (if (and (fx>=? max-depolarization 2)
                       (fx=? max-depolarization 
                             (fx+ (if (fx=? cellx basal-cellx)  2 0)
                                  (if (fx=? cellx apical-cellx) 1 0))))
                  (cons cellx predicted-cells-for-col)
                  predicted-cells-for-col))
            (if (fx<? this-colx (tm-num-columns tm))
              (each-col basal-segs apical-segs this-colx
                        (append predicted-cells-for-col predicted-cells))
              (reverse (append predicted-cells-for-col predicted-cells)))))))))
                                                                                            ;
(define (depolarize-cells                ;; TM {CellX} {CellX} Boolean ->
          tm basal-input apical-input learn)
  ;; Calculate dendrite segment activity, using the current active cells.
  (let-values ([(active-segments matching-segments potential-overlaps)
                  (compute-activity tm (list->vector basal-input) (tm-basal-pre-index tm))])
    (tm-active-basal-segments-set!    tm (sorted-segs active-segments))
    (tm-matching-basal-segments-set!  tm (sorted-segs matching-segments))
    (tm-basal-potential-overlaps-set! tm potential-overlaps))
  (let-values ([(active-segments matching-segments potential-overlaps)
                  (compute-activity tm (list->vector apical-input) (tm-apical-pre-index tm))])
    (tm-active-apical-segments-set!    tm (sorted-segs active-segments))
    (tm-matching-apical-segments-set!  tm (sorted-segs matching-segments))
    (tm-apical-potential-overlaps-set! tm potential-overlaps))
  (tm-predicted-cells-set! tm (calculate-predicted-cells tm (tm-active-basal-segments tm)
                                                            (tm-active-apical-segments tm)))
  (when learn
    (let ((set-last-used (lambda (segment)
                           (seg-last-used-set! segment (tm-iteration tm)))))
      (for-each set-last-used (tm-active-basal-segments  tm))
      (for-each set-last-used (tm-active-apical-segments tm)))
    (tm-iteration-set! tm (add1 (tm-iteration tm)))))
                                                                                            ;
(define (reset tm)                       ;; TM ->
  ;; clear history for start of new sequence of inputs
  (tm-active-cells-set! tm                  '#())
  (tm-winner-cells-set! tm                  '())
  (tm-predicted-cells-set! tm               '())
  (tm-predicted-active-cells-set! tm        '())
  (tm-active-basal-segments-set! tm         '())
  (tm-matching-basal-segments-set! tm       '())
  (tm-active-apical-segments-set! tm        '())
  (tm-matching-apical-segments-set! tm      '())
  (tm-prev-predicted-cells-set! tm          '())
  (tm-prev-apical-input-set! tm             '())
  (tm-prev-apical-growth-candidates-set! tm '())
  (tm-chosen-cell-for-column-set! tm
    (if (tm-learn-on-one-cell tm)
        (make-vector (tm-num-columns tm) #f)
        #f)))
                                                                                            ;
(define (get-active-cells tm)            ;; TM -> {CellX}
  (vector->list (tm-active-cells tm)))
                                                                                            ;
(define compute                          ;; TM {ColX} [{CellX} {CellX} [{CellX} {CellX}]] Boolean ->
  (case-lambda
  [ (tm active-columns basal-input apical-input basal-growth-candidates apical-growth-candidates learn)
      ;; general case: one time step with external, basal, and apical inputs
      (depolarize-cells tm basal-input apical-input learn)
      (tm-prev-predicted-cells-set! tm (tm-predicted-cells tm))
      (activate-cells tm active-columns basal-input apical-input 
                      basal-growth-candidates apical-growth-candidates learn) ]
                                                                                            ;
  [ (tm active-columns apical-input apical-growth-candidates learn)
      ;; Traditional TM sequence memory, with apical tiebreak
      ;; use active/winner cells from previous time step as basal input
      (tm-prev-predicted-cells-set! tm (tm-predicted-cells tm))
      (let ((prev-active-cells (vector->list (tm-active-cells tm)))
            (prev-winner-cells (tm-winner-cells tm)))
        (activate-cells tm active-columns prev-active-cells (tm-prev-apical-input tm)
                        prev-winner-cells (tm-prev-apical-growth-candidates tm) learn)
        (depolarize-cells tm (vector->list (tm-active-cells tm)) apical-input learn))
      (tm-prev-apical-input-set! tm apical-input)
      (tm-prev-apical-growth-candidates-set! tm apical-growth-candidates) ]
                                                                                            ;
  [ (tm active-columns learn)
      ;; htm-tm compatible
      (compute tm active-columns '() '() learn) ] ))
                                                                                            ;
(define (get-predicted-cells tm)         ;; TM -> {CellX}
  (tm-prev-predicted-cells tm))
                                                                                            ;
(define (get-predictive-cols tm)         ;; TM -> {ColX}
  (let loop ( (previous-colx -1) 
              (predictive-cols '()) 
              (predicted-cells (calculate-predicted-cells tm (tm-active-basal-segments tm)
                                                             (tm-active-apical-segments tm))))
    (let ((this-colx (cellxs->colx tm predicted-cells)))
      (cond [ (null? predicted-cells) (reverse predictive-cols) ]
            [ (fx=? this-colx previous-colx) 
                (loop previous-colx predictive-cols (cdr predicted-cells)) ]
            [ else (loop this-colx (cons this-colx predictive-cols) (cdr predicted-cells)) ]))))

;; === Smoke tests ===
                                                                                            ;
(define tests "")
(define failures #f)
                                                                                            ;
(define (test name)
  (set! tests (string-append tests "\n" name)))
                                                                                            ;
(define-syntax expect                    ;; ((X ... -> Y) X ...) Y -> [error]
  ;; check that function application(s) to arguments match expected values
  (lambda (x)                            
    (syntax-case x ()                                  ;; [expect ([fn args] expected ) ... ]
      [ (_ (expr expected) ...)                        ;; expr matches [fn args]
        #'(begin (let ((result expr))                  ;; eval expr just once
                   (if (equal? result expected)
                      (set! tests (string-append tests " ok"))
                      (begin 
                        (for-each display `(,tests #\newline "**" expr #\newline
                          "  expected: " ,expected #\newline 
                          "  returned: " ,result  #\newline))
                        (set! failures #t)
                        (set! tests "")))) ...)])))
                                                                                            ;
(define print-fields
  ;; Excerpt From: R. Kent Dybvig. The Scheme Programming Language, 4th Edition.
  (lambda (r)
    (unless (record? r)
      (assertion-violation 'print-fields "not a record" r))
    (let loop ([rtd (record-rtd r)])
      (let ([prtd (record-type-parent rtd)])
        (when prtd (loop prtd)))
      (let* ([v (record-type-field-names rtd)]
             [n (vector-length v)])
        (do ([i 0 (+ i 1)])
            ((= i n))
          (write (vector-ref v i))
          (display "=")
          (write ((record-accessor rtd i) r))
          (newline))))))
                                                                                            ;
(define (make-test-tm . options)
  (random-seed! 42)
  (apply make-tm 
    (append 
      (list '(32) 4)
      options
      `([activation-threshold . 3]
        [min-threshold        . 2]
        [sample-size          . 3]
        [basal-predicted-segment-decrement  . 0]
        [apical-predicted-segment-decrement . 0]
        [initial-permanence . ,(tm:permanence 0.5)] ;; so that grow-synapses creates connected synapses
      ))))
                                                                                            ;
(define (create-basal-segment tm segx)
  (create-segment tm segx (tm-basal-connections tm)))
                                                                                            ;
(define (create-basal-synapses tm seg pre-cells)
  (grow-synapses tm seg (length pre-cells) pre-cells (tm-basal-pre-index tm)))
                                                                                            ;
(test "ActivateCorrectlyPredictiveCells")
  (let* ( (tm (make-test-tm))
          (active-segment (create-basal-segment tm 4)))
    (create-basal-synapses tm active-segment '(0 1 2 3))
    (compute tm '(0) #t)
    (compute tm '(1) #t)
    [expect ([get-predicted-cells tm] '(4))
            ([get-active-cells    tm] '(4))])

(test "BurstUnpredictedColumns")
  (let* ( (tm (make-test-tm)))
    (compute tm '(0) #t)
    [expect ([get-predicted-cells tm] '())
            ([get-active-cells    tm] '(0 1 2 3))])
                                                                                            ;
(test "ZeroActiveColumns")
  (let* ( (tm (make-test-tm `[predicted-segment-decrement . ,(tm:permanence 0.02)]))
          (segment (create-basal-segment tm 4)))
    (create-basal-synapses tm segment '(0 1 2 3))
    (compute tm '(0) #t)
    [expect ((null? [get-active-cells tm]) #f)
            ((null? [tm-winner-cells tm])  #f)]
    (compute tm '() #t)
    [expect ((null? [get-active-cells tm]) #t)
            ((null? [tm-winner-cells tm])  #t)])
                                                                                            ;
(test "PredictedActiveCellsAreAlwaysWinners")
  (let* ( (tm (make-test-tm))
          (active-segment-1 (create-basal-segment tm 4))
          (active-segment-2 (create-basal-segment tm 6)))
    (create-basal-synapses tm active-segment-1 '(0 1 2))
    (create-basal-synapses tm active-segment-2 '(0 1 2))
    (compute tm '(0) '() '() #f)
    (compute tm '(1) '() '() #f)
    [expect ([tm-winner-cells tm]  '(4 6))])
                                                                                            ;
(test "ChooseOneWinnerCellInBurstingColumn")
  (let* ( (tm (make-test-tm)))
    (compute tm '(0) '() '() #f)
    [expect ((length [tm-winner-cells tm]) 1)
            ((not (member (car [tm-winner-cells tm]) '(0 1 2 3))) #f)])
                                                                                            ;
(test "ReinforceCorrectlyActiveSegments")
  (let* ( (tm (make-test-tm
                `[sample-size                 . 4]
                `[permanence-decrement        . ,(tm:permanence 0.08)]
                `[predicted-segment-decrement . ,(tm:permanence 0.02)]))
          (active-segment (create-basal-segment tm 5)))
    (create-basal-synapses tm active-segment '(0 1 2 81))
    (compute tm '(0) #t)
    (compute tm '(1) #t)
    [expect ((vector-map tm:perm [seg-synapses active-segment]) '#(6000 6000 6000 5000 4200))])
                                                                                            ;
(test "ReinforceSelectedMatchingSegmentInBurstingColumn")
  (let* ( (tm (make-test-tm
                `[permanence-decrement . ,(tm:permanence 0.08)]
                `[initial-permanence   . ,(tm:permanence 0.3)]))
          (selected-matching-segment (create-basal-segment tm 4))
          (other-matching-segment    (create-basal-segment tm 5)))
    (create-basal-synapses tm selected-matching-segment '(0 1 2 81))
    (create-basal-synapses tm other-matching-segment    '(0 1 81))
    (compute tm '(0) #t)
    (compute tm '(1) #t)
    [expect ((vector-map tm:perm [seg-synapses selected-matching-segment]) '#(4000 4000 4000 2200))])
                                                                                            ;
(test "NoChangeToNonselectedMatchingSegmentsInBurstingColumn")
  (let* ( (tm (make-test-tm
                `[permanence-decrement . ,(tm:permanence 0.08)]
                `[initial-permanence   . ,(tm:permanence 0.3)]))
          (selected-matching-segment (create-basal-segment tm 4))
          (other-matching-segment    (create-basal-segment tm 5)))
    (create-basal-synapses tm selected-matching-segment '(0 1 2 81))
    (create-basal-synapses tm other-matching-segment    '(0 1 81))
    (compute tm '(0) #t)
    (compute tm '(1) #t)
    [expect ((vector-map tm:perm [seg-synapses other-matching-segment]) '#(3000 3000 3000))])
                                                                                            ;
(test "NoChangeToMatchingSegmentsInPredictedActiveColumn")
  (let* ( (tm (make-test-tm))
          (active-segment (create-basal-segment tm 4))
          (matching-segment-on-same-cell  (create-basal-segment tm 4))
          (matching-segment-on-other-cell (create-basal-segment tm 5)))
    (create-basal-synapses tm active-segment '(0 1 2 3))
    (create-basal-synapses tm matching-segment-on-same-cell '(0 1))
    (seg-synapses-set! matching-segment-on-same-cell (vector (tm:synapse 0 3000) (tm:synapse 1 3000)))
    (create-basal-synapses tm matching-segment-on-other-cell '(0 1))
    (seg-synapses-set! matching-segment-on-other-cell (vector (tm:synapse 0 3000) (tm:synapse 1 3000)))
    (compute tm '(0) #t)
    (compute tm '(1) #t)
    [expect ([get-predicted-cells tm] '(4))
            ((vector-map tm:perm [seg-synapses matching-segment-on-same-cell])  '#(3000 3000))
            ((vector-map tm:perm [seg-synapses matching-segment-on-other-cell]) '#(3000 3000))])
                                                                                            ;
(test "NoNewSegmentIfNotEnoughWinnerCells")
  (let* ( (tm (make-test-tm
                `[sample-size                 . 2])))
    (compute tm '() #t)
    (compute tm '(0) #t)
    [expect ([tm-next-flatx tm] 0)])
                                                                                            ;
(test "NewSegmentAddSynapsesToSubsetOfWinnerCells")
  (let* ( (tm (make-test-tm
                `[sample-size        . 2]
                `[initial-permanence . ,(tm:permanence 0.21)])))
    (compute tm '(0 1 2) #t)
    (let ((prev-winner-cells [tm-winner-cells tm]))
      [expect ((length prev-winner-cells) 3)]
      (compute tm '(4) #t)
      [expect ((length [tm-winner-cells tm]) 1)]
      (let* ( (segments (vector-ref [tm-basal-connections tm] (car [tm-winner-cells tm])))
              (synapses (seg-synapses (car segments))))
        [expect ((length segments) 1)
                ((vector-length synapses) 2)]
        (vector-for-each
          (lambda (synapse)
            [expect ((tm:perm synapse) 2100)
                    ((not (member (tm:prex synapse) prev-winner-cells)) #f)])
          synapses))))
                                                                                            ;
(test "NewSegmentAddSynapsesToAllWinnerCells")
  (let* ( (tm (make-test-tm
                `[sample-size        . 4]
                `[initial-permanence . ,(tm:permanence 0.21)])))
    (compute tm '(0 1 2) #t)
    (let ((prev-winner-cells [tm-winner-cells tm]))
      [expect ((length prev-winner-cells) 3)]
      (compute tm '(4) #t)
      [expect ((length [tm-winner-cells tm]) 1)]
      (let* ( (segments (vector-ref [tm-basal-connections tm] (car [tm-winner-cells tm])))
              (synapses (seg-synapses (car segments))))
        [expect ((length segments) 1)
                ((vector-length synapses) 3)]
        (vector-for-each
          (lambda (synapse)
            [expect ((tm:perm synapse) 2100)
                    ((not (member (tm:prex synapse) prev-winner-cells)) #f)])
          synapses)
        [expect ((vector->list (vector-map tm:prex synapses)) prev-winner-cells)])))
                                                                                            ;
(test "MatchingSegmentAddSynapsesToSubsetOfWinnerCells")
  (let* ( (tm (make-tm '(32) 1
                `[basal-input-size . ,(* 32 1)]
                `[activation-threshold . 3]
                `[min-threshold . 1]
                `[sample-size . 3]
                `[predicted-segment-decrement . 0]))
          (matching-segment (create-basal-segment tm 4)))
    (create-basal-synapses tm matching-segment '(0))
    (seg-synapses-set! matching-segment (vector (tm:synapse 0 5000)))
    (compute tm '(0 1 2 3) #t)
    [expect ([tm-winner-cells tm] '(0 1 2 3))]
    (compute tm '(4) #t)
    (let* ( (synapses (seg-synapses matching-segment)))
      [expect ((vector-length synapses) 3)
              ((tm:perm (vector-ref synapses 1)) 2100)
              ((not (member (tm:prex (vector-ref synapses 1)) '(1 2 3))) #f)
              ((tm:perm (vector-ref synapses 2)) 2100)
              ((not (member (tm:prex (vector-ref synapses 2)) '(1 2 3))) #f)]))
                                                                                            ;
(test "MatchingSegmentAddSynapsesToAllWinnerCells")
  (let* ( (tm (make-tm '(32) 1
                `[basal-input-size . ,(* 32 1)]
                `[activation-threshold . 3]
                `[min-threshold . 1]
                `[sample-size . 3]
                `[predicted-segment-decrement . 0]))
          (matching-segment (create-basal-segment tm 4)))
    (create-basal-synapses tm matching-segment '(0))
    (seg-synapses-set! matching-segment (vector (tm:synapse 0 5000)))
    (compute tm '(0 1) #t)
    [expect ([tm-winner-cells tm] '(0 1))]
    (compute tm '(4) #t)
    (let* ( (synapses (seg-synapses matching-segment)))
      [expect ((vector-length synapses) 2)
              ((tm:perm (vector-ref synapses 1)) 2100)
              ((not (member (tm:prex (vector-ref synapses 1)) '(1 2 3))) #f)]))
                                                                                            ;
(test "ActiveSegmentGrowSynapsesAccordingToPotentialOverlap")
  (let* ( (tm (make-tm '(32) 1
                `[basal-input-size . ,(* 32 1)]
                `[activation-threshold . 2]
                `[min-threshold . 1]
                `[sample-size . 4]
                `[predicted-segment-decrement . 0]))
          (active-segment (create-basal-segment tm 5)))
    (create-basal-synapses tm active-segment '(0 1 2))
    (seg-synapses-set! active-segment (vector (tm:synapse 0 5000) (tm:synapse 1 5000) (tm:synapse 2 2000)))
    (compute tm '(0 1 2 3 4) #t)
    [expect ([tm-winner-cells tm] '(0 1 2 3 4))]
    (compute tm '(5) #t)
    (let* ( (synapses (seg-synapses active-segment)))
      [expect ((vector-length synapses) 4)
              ((tm:perm (vector-ref synapses 3)) 2100)
              ((not (member (tm:prex (vector-ref synapses 3)) '(3 4))) #f)]))
                                                                                            ;
(test "ActiveSegmentGrowSynapsesAccordingToPotentialOverlap")
  (let* ( (tm (make-tm '(32) 1
                `[basal-input-size . ,(* 32 1)]
                `[activation-threshold . 2]
                `[min-threshold . 1]
                `[sample-size . 4]
                `[predicted-segment-decrement . 0]))
          (active-segment (create-basal-segment tm 5)))
    (create-basal-synapses tm active-segment '(0 1 2))
    (seg-synapses-set! active-segment (vector (tm:synapse 0 5000) (tm:synapse 1 5000) (tm:synapse 2 2000)))
    (compute tm '(5) '(0 1 2 3 4) '() '(0 1 2 3 4) '() #t)
    (let* ( (synapses (seg-synapses active-segment)))
      [expect ((vector-length synapses) 4)
              ((tm:perm (vector-ref synapses 3)) 2100)
              ((not (member (tm:prex (vector-ref synapses 3)) '(3 4))) #f)]))
                                                                                            ;
  ;; flush any test failures
  (when failures
    (display tests)
    (newline))
  (flush-output-port (current-output-port))

)
