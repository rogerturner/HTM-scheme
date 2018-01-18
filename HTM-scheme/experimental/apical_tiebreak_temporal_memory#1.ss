;; === HTM-scheme Apical Tiebreak Temporal Memory Copyright 2017 Roger Turner. ===
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

  ;; Extends the HTM-scheme Temporal Memory library.
  ;; Translated from numenta htmresearch/.../apical_tiebreak_temporal_memory.py,
  ;; and htmresearch-core/.../ApicalTiebreakTemporalMemory.cpp --
  ;; see comments there for descriptions of functions and parameters.
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.

#!chezscheme

;(optimize-level 3)
                                                                                            ;
(library (apical_tiebreak_temporal_memory)
                                                                                            ;
(export
  ;; for tests in sequence_memory
  create-segment
  grow-synapses
  tm-basal-connections
  tm-basal-pre-index
  syn-perm
  syn-prex
  seg-synapses
  seg-synapses-set!
  tm-next-flatx

  tm-predicted-cells
  tm:permanence
  tm:synapse
  tm:segment
  tm:get-active-cols
  (rename
    (tm                   tm:tm)
    (make-tm              tm:construct)
    (reset                tm:reset)
    (activate-cells       tm:activate-cells)
    (depolarize-cells     tm:depolarize-cells)
    (get-active-cells     tm:get-active-cells)
    (get-winner-cells     tm:get-winner-cells)))
    
                                                                                              ;
(import 
  (except (chezscheme) add1 make-list random reset)
  (htm_prelude)
  (rename (temporal_memory)
    (tm-active-segments             tm-active-basal-segments)
    (tm-active-segments-set!        tm-active-basal-segments-set!)
    (tm-matching-segments           tm-matching-basal-segments)
    (tm-matching-segments-set!      tm-matching-basal-segments-set!)
    (tm-cells                       tm-basal-connections)
    (tm-pre-index                   tm-basal-pre-index)
    (tm-predicted-segment-decrement tm-basal-predicted-segment-decrement)))
    

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
(define-record-type tm                     ;; TM
  (parent tm:tm)
  (fields
    basal-input-size                       ;; Nat
    apical-input-size                      ;; Nat
    reduced-basal-threshold                ;; Nat
    sample-size                            ;; Nat
    ;;basal-predicted-segment-decrement    ;; Permanence [renamed tm predicted-segment-decrement]
    apical-predicted-segment-decrement     ;; Permanence
    (mutable predicted-cells)              ;; (listof CellX)
    (mutable predicted-active-cells)       ;; (listof CellX) 
    ;;(mutable active-basal-segments)      ;; (listof Segment) [renamed tm active-segments]
    ;;(mutable matching-basal-segments)    ;; (listof Segment) [renamed tm matching-segments]
    ;;basal-connections                    ;; (CellVecOf (listof Segment)) [renamed tm cells]
    (mutable basal-potential-overlaps)     ;; (vectorof Nat) [indexed by flatx]
    (mutable active-apical-segments)       ;; (listof Segment)
    (mutable matching-apical-segments)     ;; (listof Segment)
    apical-connections                     ;; Connections
    apical-pre-index                       ;; (CellVecOf {FlatX})
    (mutable apical-potential-overlaps)    ;; (vectorof Nat) [indexed by flatx]
    learn-on-one-cell                      ;; Boolean
    (mutable chosen-cell-for-column)       ;; #f|(ColVecOf CellX|#f)
    (mutable use-apical-tiebreak)          ;; Boolean
    (mutable use-apical-modulation-basal-threshold)) ;; Boolean
    
  (protocol
    (lambda (pargs->new)                   ;; Nat Nat (listof KWarg) -> TM
      (lambda (column-count cells-per-column . kwargs)
      (let ((column-dimensions (list column-count))
              (num-cells         (* column-count cells-per-column)))
          (apply (apply pargs->new (append (list column-dimensions cells-per-column) kwargs))
                 (key-word-args
                    (append `(
                              [apical-connections     . ,(make-vector num-cells '())]
                              [apical-pre-index       . ,(make-vector num-cells '())]
                              [chosen-cell-for-column . ,(if (assq 'learn-on-one-cell (car kwargs))
                                                             (make-vector column-count #f)
                                                             #f)])
                      (car kwargs))
                    attm-defaults #f)))))))
                                                                                            ;
(define attm-defaults                      ;; (listof KWarg)
  `(
    [basal-input-size                   . 0]
    [apical-input-size                  . 0]
    [reduced-basal-threshold            . 13]
    [sample-size                        . 20]
    [apical-predicted-segment-decrement . ,(tm:permanence 0.0)]
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
    [use-apical-tiebreak                . #t]
    [use-apical-modulation-basal-threshold . #t]))
    
                                                                                            ;
(define (map-segments-to-cells segs)       ;; {Segment} -> {CellX}
  (map seg-cellx segs))
                                                                                            ;
(define (cellxs->colx tm cellxs)           ;; TM {CellX} -> ColX
  ;; produce first column index from list, or terminating value if list null
  (if (null? cellxs) (tm-num-columns tm) 
      (cellx->colx tm (car cellxs))))
                                                                                            ;
(define (advance-to-cell tm cellx segs)    ;; TM CellX {Segment} -> {Segment}
  ;; produce next segment list element for given cell
  (let loop ((segs segs))
    (if (null? segs) segs
      (let ((nextseg-cellx (segs->cellx tm (cdr segs))))
        (cond [(fx=? nextseg-cellx cellx) (cdr segs)]
              [(fx<? nextseg-cellx cellx) (loop (cdr segs))]
              [ else segs])))))
                                                                                            ;
(define (advance-to-col tm colx segs)      ;; TM ColX {Segment} -> {Segment}
  ;; produce next segment list element for given column
  (let loop ((segs segs))
    (if (null? segs) segs
      (let ((nextseg-colx (cellx->colx tm (segs->cellx tm (cdr segs)))))
        (cond [(fx=? nextseg-colx colx) (cdr segs)]
              [(fx<? nextseg-colx colx) (loop (cdr segs))]
              [ else segs])))))
                                                                                            ;
(define (advance-cell-to-col               ;; TM ColX {CellX} -> {CellX}
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
(define (reset tm)                         ;; TM ->
  ;; Clear all cell and segment activity.
  (tm:reset tm)
  (tm-predicted-cells-set! tm               '())
  (tm-predicted-active-cells-set! tm        '())
  (tm-active-apical-segments-set! tm        '())
  (tm-matching-apical-segments-set! tm      '())
  (tm-basal-potential-overlaps-set! tm      '())
  (tm-apical-potential-overlaps-set! tm     '()))
                                                                                            ;

(define dc 0)

(define (depolarize-cells                  ;; TM {CellX} {CellX} Boolean ->
          tm basal-input apical-input learn)
  ;; Calculate predictions.
  (let-values ([(active-apical-segments matching-apical-segments apical-potential-overlaps)
                (calculate-apical-segment-activity 
                    tm (list->vector apical-input) (tm-apical-pre-index tm))])
    (let ((reduced-basal-threshold-cells
            (if (and learn (tm-use-apical-modulation-basal-threshold tm))
              (map-segments-to-cells active-apical-segments)
              '())))
      (let-values ([(active-basal-segments matching-basal-segments basal-potential-overlaps)
                    (calculate-basal-segment-activity tm (list->vector basal-input) 
                        (tm-basal-pre-index tm) reduced-basal-threshold-cells)])
        (tm-active-basal-segments-set!    tm (sorted-segs active-basal-segments))
        (tm-active-apical-segments-set!   tm (sorted-segs active-apical-segments))

(set! dc (add1 dc))
(when (= dc 143)
  (for-each display `( "\ndepolarize " ,dc 
      " len aas: " ,(length active-apical-segments)
      " len mas: " ,(length matching-apical-segments)
      " len aps: " ,(vector-length apical-potential-overlaps)
      " len abs: " ,(length active-basal-segments)
      " len mbs: " ,(length matching-basal-segments)
      " len bps: " ,(vector-length basal-potential-overlaps)
      ))
  )
  
        (tm-predicted-cells-set! tm 
          (calculate-predicted-cells tm (tm-active-basal-segments tm) (tm-active-apical-segments tm)))
        (tm-matching-basal-segments-set!  tm (sorted-segs matching-basal-segments))
        (tm-matching-apical-segments-set! tm (sorted-segs matching-apical-segments))
        (tm-basal-potential-overlaps-set! tm basal-potential-overlaps)
        (tm-apical-potential-overlaps-set! tm apical-potential-overlaps)))))
                                                                                          ;
(define (activate-cells                    ;; TM {ColX} {CellX} {CellX} {CellX} {CellX} Boolean ->
          tm active-columns basal-reinforce-candidates apical-reinforce-candidates
          basal-growth-candidates apical-growth-candidates learn)
  ;; Activate cells in the specified columns, using the result of the previous
  ;; 'depolarizeCells' as predictions. Then learn.
  ;; Implemented by merging active columns, predicted cells, and matching segments,
  ;; building active/winner/predicted-active cell lists.
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
          
          #;(for-each display `( "activate-cells:" #\newline
            "active " ,(take 10 (reverse active-cells)) #\newline
            "winner " ,(take 10 (reverse winner-cells)) #\newline
            "predac " ,(take 10 predicted-active-cells) #\newline
            ))
            
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
(define (calculate-apical-segment-activity ;; TM (vectorof CellX) (CellVecOf {FlatX}) -> {Segment} {Segment} (vectorof Nat)
          tm active-input pre-index)
  ;; Calculate the active and matching apical segments for this timestep.
  (compute-activity tm active-input pre-index 0 '()))
                                                                                            ;
(define (calculate-basal-segment-activity  ;; TM (vectorof CellX) (CellVecOf {FlatX}) {CellX} -> {Segment} {Segment} (vectorof Nat)
          tm active-input pre-index reduced-basal-threshold-cells)
  ;; Calculate the active and matching basal segments for this timestep.
  (compute-activity tm active-input pre-index 
                    (tm-reduced-basal-threshold tm) reduced-basal-threshold-cells))
                                                                                            ;
#;(define (calculate-predicted-cells         ;; TM {Segment} {Segment} -> {CellX}
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

(define (list-intersect l1 l2)
  (bitwise->list (bitwise-and (list->bitwise l1) (list->bitwise l2))))
                                                                                            ;
(define (list-difference l1 l2)            ;; {Nat} {Nat} -> {Nat}
  (bitwise->list (bitwise-and (list->bitwise l1) (bitwise-not (list->bitwise l2)))))
                                                                                            ;
(define (cols->cells l n)                  ;; (listof ColX) Nat -> (listof CellX)
  (apply append (map
                  (lambda (e)
                    (build-list n (lambda (i) (+ (* e n) i))))
                  l)))

(define ncpc 0)

(define (dl l)
  (take 10 l))
  
(define (calculate-predicted-cells         ;; TM {Segment} {Segment} -> {CellX}
          tm active-basal-segments active-apical-segments)
  ;; Calculate the predicted cells, given the set of active segments.
  (let* ( (cells-for-basal-segments  (map-segments-to-cells active-basal-segments))
          (cells-for-apical-segments (map-segments-to-cells active-apical-segments))
          (fully-depolarized-cells         ;; cells with both types of segments active
            (list-intersect  cells-for-basal-segments cells-for-apical-segments))
          (partly-depolarized-cells        ;; cells with basal only
            (list-difference cells-for-basal-segments fully-depolarized-cells))
          (cpc  (tm-cells-per-column tm))
          (/cpc (lambda (cellx) (quotient cellx cpc)))
          (inhibited-columns               ;; columns with both
            (map /cpc fully-depolarized-cells))
          (inhibited-cells (cols->cells inhibited-columns cpc))
          (inhibited-mask  (list->bitwise inhibited-cells))
          (predicted-cells
            (list-sort <
              (append fully-depolarized-cells
                      (bitwise->list (bitwise-and (bitwise-not inhibited-mask)
                                                  (list->bitwise partly-depolarized-cells)))))))
  (set! ncpc (add1 ncpc))
  (when (= ncpc 143)
    (for-each display `( "\nc-p-c " ,ncpc
      "\ncbs: " ,(dl cells-for-basal-segments)
      "\ncas: " ,(dl cells-for-apical-segments)
      "\nfdc: " ,(dl fully-depolarized-cells)
      "\npdc: " ,(dl partly-depolarized-cells)
      "\nfdl: " ,(dl inhibited-columns)))
    (display "\ni-c: ")
    (do ((i 0 (+ 32 i))) ((> i 256))
      (for-each display `( 
        ,(list-ref inhibited-cells i) "-" ,(+ 31 (list-ref inhibited-cells i)) " ")))
    (for-each display `(
      "\np-c: " ,(dl predicted-cells))))

    predicted-cells

    #;(list-sort <
      (append fully-depolarized-cells
              (bitwise->list (bitwise-and (bitwise-not inhibitedMask)
                                          (list->bitwise partly-depolarized-cells)))))))
                                                                                            ;
#;(define (calculate-predicted-cells         ;; TM {Segment} {Segment} -> {CellX}
          tm active-basal-segments active-apical-segments)
  ;; Calculate the predicted cells, given the set of active segments.
  (let* ( (cells-for-basal-segments  (map-segments-to-cells active-basal-segments))
          (cells-for-apical-segments (map-segments-to-cells active-apical-segments))
          (fully-depolarized-cells         ;; cells with both types of segments active
            (list-intersect  cells-for-basal-segments cells-for-apical-segments))
          (partly-depolarized-cells        ;; cells with basal only
            (list-difference cells-for-basal-segments fully-depolarized-cells))
          (cpc  (tm-cells-per-column tm))
          (/cpc (lambda (cellx) (quotient cellx cpc)))
          (inhibited-list                  ;; columns 
            (list-intersect (map /cpc partly-depolarized-cells) 
                            (map /cpc fully-depolarized-cells)))
          (inhibited-mask (list->bitwise (list-multiples inhibited-list cpc)))
          (predicted-cells
            (list-sort <
              (append fully-depolarized-cells
                      (bitwise->list (bitwise-and (bitwise-not inhibited-mask)
                                                  (list->bitwise partly-depolarized-cells)))))))
    predicted-cells))
                                                                                            ;
(define (learn-on-cell                     ;; TM Connections CellX {Seg} {Seg} ... ->
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
(define (activate-predicted-column         ;; TM {CellX} {CellX} {CellX} ColX ... -> {CellX} {CellX} {CellX}
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
(define (segment-with-max-overlaps         ;; TM CellX {Segment} (vectorof Nat) -> {Segment}
          tm colx segments overlaps)
  (let loop ((segments segments) (best-seg '()) (max-overlap (least-fixnum)))
    (if (fx=? colx (segs->colx tm segments))
        (let ((overlap (vector-ref overlaps (seg-flatx (car segments)))))
          (if (fx>? overlap max-overlap)
              (loop (cdr segments) segments overlap)
              (loop (cdr segments) best-seg max-overlap)))
        best-seg)))
                                                                                            ;
(define (burst-column                      ;; TM {CellX} {CellX} ColX ... -> {CellX} {CellX}
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
(define (get-active-cells tm)              ;; TM -> {CellX}
  (vector->list (tm-active-cells tm)))
                                                                                            ;
(define (get-predicted-active-cells tm)    ;; TM -> {CellX}
  (tm-predicted-active-cells tm))
                                                                                            ;
(define (get-winner-cells tm)              ;; TM -> {CellX}
  (tm-winner-cells tm))
                                                                                            ;
)