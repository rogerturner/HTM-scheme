;; ===  HTM-scheme Extended Temporal Memory Copyright 2017 Roger Turner.  ===
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Based on code from Numenta Platform for Intelligent Computing (NuPIC) ;;
  ;; which is Copyright (C) 2014-2016, Numenta, Inc.                       ;;
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

  ;; Extends the HTM-scheme Temporal Memory library htm-tm.ss -- extends TM record
  ;; and uses htm-tm procedures where applicable.
  ;; Translated from NuPIC ExtendedTemporalMemory.cpp, see there for more info.
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.

#!chezscheme

(optimize-level 3)
                                                                                            ;
(library (libraries htm-etm)
                                                                                            ;
(export
  tm:initialize
  tm:permanence
  tm:get-active-cols
  (rename
    (make-tm  tm:constructor)
    (compute  tm:compute)
    (reset    tm:reset)
    (get-predictive-cols tm:get-predictive-cols)
    (get-predicted-cells tm:get-predicted-cells)))
                                                                                            ;
(import 
  (except (chezscheme) add1 make-list random reset)
  (libraries htm-prelude)
  (rename (libraries htm-tm)
    (tm-active-segments        tm-active-basal-segments)
    (tm-active-segments-set!   tm-active-basal-segments-set!)
    (tm-matching-segments      tm-matching-basal-segments)
    (tm-matching-segments-set! tm-matching-basal-segments-set!)
    (tm-cells                  tm-basal-connections)
    (tm-pre-index              tm-basal-pre-index)
    ))

;; === Temporal Memory Types ===
                                                                                            ;
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
    (mutable predicted-cells)            ;; (listof CellX)
    (mutable predicted-active-cells)     ;; (listof CellX) 
    ;;(mutable active-basal-segments)    ;; (listof Segment) renamed tm active-segments
    ;;(mutable matching-basal-segments)  ;; (listof Segment) renamed tm matching-segments
    ;;(basal-connections)                ;; (CellVecOf (listof Segment)) renamed tm cells
    (mutable basal-overlaps)             ;; (vectorof ) [flatx]
    (mutable basal-potential-overlaps)   ;; (vectorof ) [flatx]
    (mutable active-apical-segments)     ;; (listof Segment)
    (mutable matching-apical-segments)   ;; (listof Segment)
    apical-connections                   ;; Connections
    apical-pre-index
    (mutable apical-overlaps)            ;; (vectorof ) [flatx]
    (mutable apical-potential-overlaps)  ;; (vectorof ) [flatx]
    learn-on-one-cell                    ;; Boolean
    (mutable chosen-cell-for-column)     ;; (ColVecOf CellX)
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
                      [chosen-cell-for-column . ,(make-vector num-columns #f)])
                      args)
                    tm-defaults)))))))
                                                                                            ;
(define tm-defaults                      ;; (listof KWarg)
  `(
    [basal-input-size            . 0]
    [apical-input-size           . 0]
    [sample-size                 . 100]
    [predicted-cells             . ()]
    [predicted-active-cells      . ()]
    [basal-overlaps              . #()]
    [basal-potential-overlaps    . #()]
    [active-apical-segments      . ()]
    [matching-apical-segments    . ()]
    [apical-connections          . #()]
    [apical-pre-index            . #()]
    [apical-overlaps             . #()]
    [apical-potential-overlaps   . #()]
    [learn-on-one-cell           . #f]
    [chosen-cell-for-column      . #()]
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

;; === Extended Temporal Memory Algorithm ===
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
        (if (fx=? (cellx->colx tm cellx) colx)
          (begin
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
          (let* ( (predicted-cells (cdr predicted-cells))
                  (next-cellx (if (null? predicted-cells) (tm-num-cells tm)
                                                          (car predicted-cells))))
            (let ((active-basal-segs    (advance-to-cell tm next-cellx active-basal-segs))
                  (matching-basal-segs  (advance-to-cell tm next-cellx matching-basal-segs))
                  (active-apical-segs   (advance-to-cell tm next-cellx active-apical-segs))
                  (matching-apical-segs (advance-to-cell tm next-cellx matching-apical-segs)))
              (loop predicted-cells active-basal-segs matching-basal-segs active-apical-segs matching-apical-segs
                    active-cells winner-cells predicted-active-cells))))))))
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
          apical-input apical-potential-overlaps apical-growth-candidates
          learn)
  ;; choose winner cell and learn on it; add all cells in column to active cells
  (let* ( (last (fx- (fx* (add1 colx) (tm-cells-per-column tm)) 1))
          (cells-for-column (build-list (tm-cells-per-column tm)
                                        (lambda (x) (fx- last x))))
          (chosen-cell (vector-ref (tm-chosen-cell-for-column tm) colx)))
    (let-values ([(basal-candidates winner-cell)
      (if (and (tm-learn-on-one-cell tm) chosen-cell)
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
(define (activate-cells                  ;; TM {ColX} Boolean ->
          tm active-columns basal-input apical-input
          basal-growth-candidates apical-growth-candidates learn)
  ;; step through columns in merge of active columns, predicted cells, and matching 
  ;; basal segments, building active/winner/predicted-active cell lists; this blends
  ;; in the required aspects of iterGroupBy logic from GroupBy.hpp
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
          (matbsegs-colx  (segs->colx   tm matching-basal-segs)))
      (let ((next (lambda (basal? actcells wincells predactcells)
                    (let* ( (active-cols 
                              (if (or basal? (null? active-cols)) active-cols 
                                  (cdr active-cols)))
                            (matching-basal-segs
                              (if basal?
                                  (skip-col tm matbsegs-colx matching-basal-segs)
                                  (advance-to-col tm (colxs->colx tm active-cols) matching-basal-segs)))
                            (next-colx (fxmin (colxs->colx tm active-cols)
                                              (segs->colx  tm matching-basal-segs))))
                      (loop active-cols
                            (advance-cell-to-col tm next-colx predicted-cells) 
                            (advance-to-col tm next-colx active-basal-segs) 
                            matching-basal-segs 
                            (advance-to-col tm next-colx active-apical-segs) 
                            (advance-to-col tm next-colx matching-apical-segs) 
                            actcells wincells predactcells))))
            (colx (fxmin actcols-colx predcells-colx matbsegs-colx)))
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
                        basal-input basal-growth-candidates apical-input apical-growth-candidates learn)])
                    (next #f active-cells winner-cells predicted-active-cells))
                  ;; else: active but not predicted
                  (let-values 
                    ([(active-cells winner-cells)
                      (burst-column tm active-cells winner-cells colx
                        active-basal-segs matching-basal-segs active-apical-segs matching-apical-segs
                        basal-input (tm-basal-potential-overlaps tm) basal-growth-candidates 
                        apical-input (tm-apical-potential-overlaps tm) apical-growth-candidates learn)])
                    (next #f active-cells winner-cells predicted-active-cells)))
            ;; else: not active but colx in predicted|*-segs
            (begin
              (when (and learn (not (null? matching-basal-segs)))
                  (punish-predicted-column tm 
                    (segs->colx tm matching-basal-segs) matching-basal-segs (list->vector basal-input)))
              (next #t active-cells winner-cells predicted-active-cells))))))))
                                                                                            ;
(define (depolarize-cells                ;; TM {CellX} {CellX} Bool ->
          tm basal-input apical-input learn)
  ;; Calculate dendrite segment activity, using the current active cells.
  ;; Set predicted-cells
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
  ;; if MIN_PREDICTIVE_THRESHOLD is 2 then predicted cells are those with active basal
  ;; segments but no apical segments
  (let loop ( (basal-segs  (tm-active-basal-segments  tm))
              (apical-segs (tm-active-apical-segments tm))
              (predicted-cells '()))
    (if (null? basal-segs)
      (tm-predicted-cells-set! tm (reverse predicted-cells))
      (let* ( (cellx (segs->cellx tm basal-segs))
              (apical-segs (advance-to-cell tm cellx apical-segs)))
        (if (fx=? cellx (segs->cellx tm apical-segs))
            (loop (advance-to-cell tm (add1 cellx) basal-segs) apical-segs predicted-cells)
            (loop (cdr basal-segs) apical-segs (cons cellx predicted-cells))))))
  (when learn
    (let ((set-last-used (lambda (segment)
                           (seg-last-used-set! segment (tm-iteration tm)))))
      (for-each set-last-used (tm-active-basal-segments  tm))
      (for-each set-last-used (tm-active-apical-segments tm)))
    (tm-iteration-set! tm (add1 (tm-iteration tm)))))
                                                                                            ;
(define compute                          ;; TM {ColX} [{CellX} {CellX} [{CellX} {CellX}]] Bool ->
  (case-lambda
  [(tm active-columns basal-input apical-input basal-growth-candidates apical-growth-candidates learn)
    ;; general case: one time step with external, basal, and apical inputs
    (depolarize-cells tm basal-input apical-input learn)
    (activate-cells tm (list-sort fx<? active-columns) basal-input apical-input 
                    basal-growth-candidates apical-growth-candidates learn)]
  [(tm active-columns apical-input apical-growth-candidates learn)
    ;; use active/winner cells from previous time step as basal input
    (compute tm active-columns (vector->list (tm-active-cells tm)) apical-input
           (tm-winner-cells tm) apical-growth-candidates learn)]
  [(tm active-columns learn)
    ;; htm-tm compatible
    (compute tm active-columns '() '() learn)]))
                                                                                            ;
(define (reset tm)                       ;; TM ->
  ;; clear history for start of new sequence of inputs
  (tm-active-cells-set! tm             '#())
  (tm-winner-cells-set! tm             '())
  (tm-predicted-cells-set! tm          '())
  (tm-predicted-active-cells-set! tm   '())
  (tm-active-basal-segments-set! tm    '())
  (tm-matching-basal-segments-set! tm  '())
  (tm-active-apical-segments-set! tm   '())
  (tm-matching-apical-segments-set! tm '())
  (tm-chosen-cell-for-column-set! tm 
    (make-vector (tm-num-columns tm) #f)))
                                                                                            ;
(define (get-active-cells tm)            ;; TM -> {CellX}
  (vector->list (tm-active-cells tm)))
                                                                                            ;
(define (get-predicted-cells tm)         ;; TM -> {CellX}
  (tm-predicted-cells tm))
                                                                                            ;
(define (get-predictive-cols tm)         ;; TM -> (listof ColX)
  ;(depolarize-cells tm (vector->list (tm-active-cells tm)) '() #f)  ;; set predicted-cells
  (let loop ((previous-colx -1) (predictive-cols '()) (predicted-cells (tm-predicted-cells tm)))
    (let ((this-colx (cellxs->colx tm predicted-cells)))
      (cond [ (null? predicted-cells) (reverse predictive-cols) ]
            [ (fx=? this-colx previous-colx) 
                (loop previous-colx predictive-cols (cdr predicted-cells)) ]
            [ else (loop this-colx (cons this-colx predictive-cols) (cdr predicted-cells)) ]))))

;; === Smoke tests ===
                                                                                            ;
(define (tm:initialize)
  ;; flush any test failures
  (flush-output-port (current-output-port)))
                                                                                            ;
(define-syntax expect                    ;; ((X ... -> Y) X ...) Y -> [error]
  ;; check that function application(s) to arguments match expected values
  (lambda (x)                            
    (syntax-case x ()                                  ;; [expect ([fn args] expected ) ... ]
      [ (_ (expr expected) ...)                        ;; expr matches [fn args]
        #'(begin (let ((result expr))                  ;; eval expr just once
                   (unless (equal? result expected)
                       (for-each display `(#\newline "**" expr #\newline
                         "  expected: " ,expected #\newline 
                         "  returned: " ,result  #\newline)))) ...)])))
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
      `([basal-input-size . ,(* 32 4)]
        [activation-threshold . 3]
        [min-threshold . 2]
        [sample-size . 3]
        [predicted-segment-decrement . 0]
        [initial-permanence . ,(tm:permanence 0.5)] ;; so that grow-synapses creates connected synapses
      ))))
                                                                                            ;
(define (create-basal-segment tm segx)
  (create-segment tm segx (tm-basal-connections tm)))
                                                                                            ;
(define (create-basal-synapses tm seg pre-cells)
  (grow-synapses tm seg (length pre-cells) pre-cells (tm-basal-pre-index tm)))
                                                                                            ;
;; ActivateCorrectlyPredictiveCells
  (let* ( (tm (make-test-tm))
          (active-segment (create-basal-segment tm 4)))
    (create-basal-synapses tm active-segment '(0 1 2 3))
    (compute tm '(0) #t)
    (compute tm '(1) #t)
    [expect ([get-predicted-cells tm] '(4))
            ([get-active-cells    tm] '(4))])
                                                                                            ;
;; ActivateCorrectlyPredictiveCells
  (let* ( (tm (make-test-tm))
          (active-segment-1 (create-basal-segment tm 4))
          (active-segment-2 (create-basal-segment tm 8)))
    (create-basal-synapses tm active-segment-1 '(0 1 2 3))
    (create-basal-synapses tm active-segment-2 '(0 1 2 3))
    (compute tm '(0) #t)
    [expect ((let ((predicted [get-predictive-cols tm]))
                (or (equal? predicted '(1 2))
                    (equal? predicted '())))   #t )]
    (compute tm '(1 2) #t)
    [expect ([get-active-cells    tm] '(4 8))])
                                                                                            ;
;; BurstUnpredictedColumns
  (let* ( (tm (make-test-tm)))
    (compute tm '(0) #t)
    [expect ([get-predicted-cells tm] '())
            ([get-active-cells    tm] '(0 1 2 3))])
                                                                                            ;
;; ZeroActiveColumns
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
;; PredictedActiveCellsAreAlwaysWinners
  (let* ( (tm (make-test-tm))
          (active-segment-1 (create-basal-segment tm 4))
          (active-segment-2 (create-basal-segment tm 6)))
    (create-basal-synapses tm active-segment-1 '(0 1 2))
    (create-basal-synapses tm active-segment-2 '(0 1 2))
    (compute tm '(0) '() '() #f)
    (compute tm '(1) '() '() #f)
    [expect ([tm-winner-cells tm]  '(4 6))])
                                                                                            ;
;; ChooseOneWinnerCellInBurstingColumn
  (let* ( (tm (make-test-tm)))
    (compute tm '(0) '() '() #f)
    [expect ((length [tm-winner-cells tm]) 1)
            ((not (member (car [tm-winner-cells tm]) '(0 1 2 3))) #f)])
                                                                                            ;
;; ReinforceCorrectlyActiveSegments
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
;; ReinforceSelectedMatchingSegmentInBurstingColumn
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
;; NoChangeToNonselectedMatchingSegmentsInBurstingColumn
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
;; NoChangeToMatchingSegmentsInPredictedActiveColumn
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
;; NoNewSegmentIfNotEnoughWinnerCells
  (let* ( (tm (make-test-tm
                `[sample-size                 . 2])))
    (compute tm '() #t)
    (compute tm '(0) #t)
    [expect ([tm-next-flatx tm] 0)])
                                                                                            ;
;; NewSegmentAddSynapsesToSubsetOfWinnerCells
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
;; NewSegmentAddSynapsesToAllWinnerCells
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
;; MatchingSegmentAddSynapsesToSubsetOfWinnerCells
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
;; MatchingSegmentAddSynapsesToAllWinnerCells
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
;; ActiveSegmentGrowSynapsesAccordingToPotentialOverlap
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
    
  (newline)
  
  )
