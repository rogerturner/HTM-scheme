#!chezscheme

;; ====== HTM-scheme Spatial Pooler Copyright 2017-2020 Roger Turner. ======
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

  ;; Translated from NuPIC spatial_pooler.py, see comments there for more info.
  ;; Max 2 dimensions, wraparound for 1D only, no parameter consistency checking.
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.

(library (HTM-scheme HTM-scheme algorithms spatial_pooler)
                                                                                            ;
(export
  sp
  make-sp
  compute
  sp-num-columns
  ;; following exports enable override of compute/adapt-synapses
  sp-iteration-num
  sp-iteration-num-set!
  sp-iteration-learn-num
  sp-iteration-learn-num-set!
  calculate-overlap
  sp-connected-counts
  inhibit-columns
  update-duty-cycles
  bump-up-weak-columns
  update-boost-factors
  sp-update-period
  update-inhibition-radius
  sp-inhibition-radius-set!
  update-min-duty-cycles
  sp-stimulus-threshold
  raise-permanence-to-threshold
  increase-perm
  decrease-perm
  sp-syn-perm-connected
  sp-syn-perm-inactive-dec
  sp-syn-perm-active-inc
  sp-syn-perm-trim-threshold
  sp-boost-factors
  for-each-segment
  )
                                                                                            ;
(import 
  (except (chezscheme) add1 make-list random reset)
  (HTM-scheme HTM-scheme algorithms htm_prelude)
  (HTM-scheme HTM-scheme algorithms htm_concept))
    
;; -- Spatial Pooler Types --
                                                                                            ;
;; InputVec     = Number (Bignum), interpreted bitwise as (vectorof Boolean)
;; InputX       = Nat, index of bit in InputVec
;; Overlap      = Flonum (to allow normalisation/fatigue of overlap counts)
;; DutyCycle    = Flonum
;; BoostFactor  = Flonum
;; ColumnX      = Nat, column index of Segment/DutyCycle/...
;; (ColVecOf X) = Vector of Xs indexed by ColumnX
;; Columns      = (ColVecOf Segment), (ColVecOf DutyCycle) ... the columns, part of SP
;; SP           = Record, sp parameters accessed as (sp-param-name sp)
                                                                                            ;
#| "_potentialPools" and "_permanences" are combined in segments (vectors of synapses)
   connected-synapses are a bitwise mask of synapses with perm >= connected |#

;; -- Spatial Pooler data --
                                                                                            ;
(define-record-type sp                   ;; SP
  (fields
    input-dimensions                          ;; (listof Nat)
    column-dimensions                         ;; (listof Nat)
    (mutable potential-radius)                ;; Nat
    (mutable potential-pct)                   ;; Number
    (mutable global-inhibition)               ;; Boolean
    (mutable num-active-columns-per-inh-area) ;; Nat
    (mutable local-area-density)              ;; Number
    (mutable stimulus-threshold)              ;; Nat
    (mutable inhibition-radius)               ;; Nat
    (mutable duty-cycle-period)               ;; Nat
    boost-strength                            ;; Number
    wrap-around                               ;; Boolean
    (mutable iteration-num)                   ;; Nat
    (mutable iteration-learn-num)             ;; Nat
    (mutable update-period)                   ;; Nat
    (mutable syn-perm-trim-threshold)         ;; Permanence
    (mutable syn-perm-active-inc)             ;; Permanence
    (mutable syn-perm-inactive-dec)           ;; Permanence
    (mutable syn-perm-below-stimulus-inc)     ;; Permanence
    (mutable syn-perm-connected)              ;; Permanence
    (mutable min-pct-overlap-duty-cycles)     ;; DutyCycle
    (mutable boost-factors)                   ;; (colVecOf BoostFactor)
    (mutable overlap-duty-cycles)             ;; (ColVecOf DutyCycle)
    (mutable active-duty-cycles)              ;; (ColVecOf DutyCycle)
    (mutable min-overlap-duty-cycles)         ;; (ColVecOf DutyCycle)
    (mutable potential-pools)                 ;; (ColVecOf Segment)
    (mutable connected-synapses)              ;; (ColVecOf InputVec)
    (mutable connected-counts))               ;; (ColVecOf Nat)
  (protocol
    (lambda (new)                        ;; (listof KWarg) -> SP
      (lambda (kwargs)
        (let* (
            [sp          (apply new (key-word-args kwargs sp-defaults))]
            [num-columns (sp-num-columns sp)])
          (sp-syn-perm-below-stimulus-inc-set! sp (fxdiv (sp-syn-perm-connected sp) 10))
          (sp-syn-perm-trim-threshold-set!     sp (fxmax (fxdiv (sp-syn-perm-active-inc sp) 2)
                                                         (sp-syn-perm-below-stimulus-inc sp)))
          (sp-potential-pools-set!             sp (build-vector num-columns (lambda (cx)
                                                      (init-permanence sp (map-potential sp cx)))))
          (sp-connected-synapses-set!          sp (build-vector num-columns (lambda (cx)
                                                      (connected-inputs sp cx))))
          (sp-connected-counts-set!            sp (make-vector (sp-num-columns sp)))
          (for-each-segment sp (lambda (segment cx)
              (raise-permanence-to-threshold sp segment))
            (build-list num-columns id))
          (sp-boost-factors-set!               sp (make-vector num-columns 1.0))
          (sp-overlap-duty-cycles-set!         sp (make-vector num-columns 0.0))
          (sp-active-duty-cycles-set!          sp (make-vector num-columns 0.0))
          (sp-min-overlap-duty-cycles-set!     sp (make-vector num-columns 0.0))
          (sp-inhibition-radius-set!           sp (update-inhibition-radius sp))
          sp)))))

(define sp-defaults `(                   ;; (listof [Keyword . X])
    (input-dimensions                . (32 32))
    (column-dimensions               . (64 64))
    (potential-radius                . 16)
    (potential-pct                   . 0.5)
    (global-inhibition               . #f)
    (num-active-columns-per-inh-area . 10)
    (local-area-density              . -1.0)
    (stimulus-threshold              . 0)
    (inhibition-radius               . 0)
    (duty-cycle-period               . 1000)
    (boost-strength                  . 0.0)
    (wrap-around                     . #t)
    (iteration-num                   . 0)
    (iteration-learn-num             . 0)
    (update-period                   . 50)
    (syn-perm-trim-threshold         . ,(perm  0.025))
    (syn-perm-active-inc             . ,(perm  0.050))
    (syn-perm-inactive-dec           . ,(perm  0.008))
    (syn-perm-below-stimulus-inc     . ,(perm  0.010))
    (syn-perm-connected              . ,(perm  0.100))
    (min-pct-overlap-duty-cycles     . 0.001)
    (boost-factors                   . #())
    (overlap-duty-cycles             . #())
    (active-duty-cycles              . #())
    (min-overlap-duty-cycles         . #())
    (potential-pools                 . #())
    (connected-synapses              . #())
    (connected-counts                . #())))
                                                                                            ;
(define (sp-num-columns sp)              ;; SP -> Nat
  (apply fx* (sp-column-dimensions sp)))
                                                                                            ;
(define (sp-num-inputs  sp)              ;; SP -> Nat
  (apply fx* (sp-input-dimensions  sp)))
                                                                                            ;
(define init-connected-pct 0.5)
                                                                                            ;
;; -- Compute --
                                                                                            ;
(define (compute sp input-vector learn)  ;; SP InputVec Boolean -> (listof ColumnX)
  ;; produce active columns from input; optionally update sp if learning
  (sp-iteration-num-set! sp (fx1+ (sp-iteration-num sp)))
  (when learn (sp-iteration-learn-num-set! sp (fx1+ (sp-iteration-learn-num sp))))
  (let* (
      [overlaps         (calculate-overlap sp input-vector)]
      [boosted-overlaps (if learn
                          (vector-map fl* overlaps (sp-boost-factors sp))
                          overlaps)]
      [active-columns   (inhibit-columns sp boosted-overlaps)])
    (when learn
      (adapt-synapses       sp input-vector active-columns)
      (update-duty-cycles   sp overlaps active-columns)
      (bump-up-weak-columns sp)
      (update-boost-factors sp)
      (when (fxzero? (fxmod (sp-iteration-num sp) (sp-update-period sp)))
        (sp-inhibition-radius-set! sp (update-inhibition-radius sp))
        (update-min-duty-cycles sp)))
    active-columns))
                                                                                            ;
;; -- Duty cycles --
                                                                                            ;
(define (update-min-duty-cycles sp)      ;; SP ->
  ;; recompute minimum overlap duty cycles as percentage of maximum
  (if (or (sp-global-inhibition sp) (fx>? (sp-inhibition-radius sp) (sp-num-inputs sp)))
      (update-min-duty-cycles-global sp))
      (update-min-duty-cycles-local  sp))
                                                                                              ;
(define (flvector-max vec)               ;; (Vectorof Flonum) -> Flonum
  (vector-fold-left flmax (vector-ref vec 0) vec))
                                                                                              ;
(define (update-min-duty-cycles-global sp);; SP ->
  ;; produce minodc as percent of max odc of all columns
  (vector-fill! (sp-min-overlap-duty-cycles sp)
    (fl* (sp-min-pct-overlap-duty-cycles sp)
       (flvector-max (sp-overlap-duty-cycles sp)))))
                                                                                            ;
(define (update-min-duty-cycles-local sp);; SP ->
  ;; produce minodc of column as percent of max odc of columns in its neighborhood
  (for-each-column sp (lambda (cx)
    (let* (
        [neighborhood     
          (neighborhood cx (sp-inhibition-radius sp) (sp-column-dimensions sp) (sp-wrap-around sp))]
        [max-overlap-duty (flvector-max (vector-refs (sp-overlap-duty-cycles sp) neighborhood))])
      (vector-set! (sp-min-overlap-duty-cycles sp) cx
        (fl* max-overlap-duty (sp-min-pct-overlap-duty-cycles sp)))))))
                                                                                            ;
(define (update-duty-cycles              ;; SP (ColVecOf Overlap) (listof ColumnX) ->
          sp overlaps active-columns)
  ;; recompute duty cycles from this input's overlaps and resulting active columns
  (let ((overlap-array (vector-map (lambda (ov) (if (positive? ov) 1.0 0.0))
                                    overlaps))
        (active-array  (make-vector (sp-num-columns sp) 0.0)))
    (for-each 
      (lambda (cx) (vector-set! active-array cx 1.0)) 
      active-columns)
    (sp-overlap-duty-cycles-set! sp 
      (update-duty-cycles-helper (sp-overlap-duty-cycles sp) overlap-array (sp-duty-cycle-period sp)))
    (sp-active-duty-cycles-set!  sp 
      (update-duty-cycles-helper (sp-active-duty-cycles  sp) active-array  (sp-duty-cycle-period sp)))))
                                                                                            ;
(define (update-duty-cycles-helper       ;; (ColVecOf DutyCycle) (ColVecOf Number) Nat -> (ColVecOf DutyCycle)
          duty-cycles new-input period)  
  ;; produce updated duty cycle estimate
  (vector-map (lambda (dc ni)
                (/ (+ (* dc (- period 1)) ni) period))
    duty-cycles new-input))
                                                                                            ;
;; -- Inhibition radius --
                                                                                            ;
(define (update-inhibition-radius sp)    ;; SP -> Nat
  ;; produce new inhibition radius
  (if (sp-global-inhibition sp)
      (apply fxmax (sp-column-dimensions sp))       ;; global
      (let* ( (avg-connected-span 
                (vector-average                     ;; local
                  (build-vector (sp-num-columns sp)
                                (lambda (cx)
                                  (if (null? (cdr (sp-column-dimensions sp)))
                                    (avg-connected-span-1d sp cx)
                                    (avg-connected-span-2d sp cx))))))
              (diameter (* avg-connected-span (avg-columns-per-input sp))))
        (int<- (max 1.0 (/ (- diameter 1) 2))))))
                                                                                            ;
(define (avg-columns-per-input sp)       ;; SP -> Number
  ;; produce mean over dimensions of num-cols / num-inputs
  (list-average (map / (sp-column-dimensions sp) (sp-input-dimensions sp))))
                                                                                            ;
(define (avg-connected-span-1d sp cx)    ;; SP ColumnX -> Nat
  ;; produce width of range of input indices of connected synapses for column
  (let ((connected (vector-ref (sp-connected-synapses sp) cx)))
    (if (zero? connected)
        0
        (bitwise-span connected))))
                                                                                            ;
(define (avg-connected-span-2d sp cx)    ;; SP ColumnX -> Number
  ;; produce average width of range of input indices of connected synapses for column
  (let ((connected (vector-ref (sp-connected-synapses sp) cx))
        (nrows     (car  (sp-input-dimensions sp)))
        (ncols     (cadr (sp-input-dimensions sp))))
    (let next-row ((rx 0) (l '(0)))
      (cond [(fx=? rx nrows) (/ (fx+ (apply fxmax l) (fx- (length l) 1)) 2.0)]
            [else
              (let* ((rowstart (fx* rx ncols))
                     (row (bitwise-bit-field connected rowstart (fx+ rowstart ncols))))
                (next-row (fx1+ rx) 
                          (if (zero? row) l
                            (cons (bitwise-span row) l))))]))))
                                                                                            ;
;; -- Adapt synapses --
                                                                                            ;
(define (increase-perm synapse by)       ;; Synapse Permanence -> Synapse
  ;; produce synapse with incremented permanence (clipped to max)
  (make-syn (syn-prex synapse) (clip-max (fx+ (syn-perm synapse) by))))
                                                                                            ;
(define (decrease-perm synapse by trim)  ;; Synapse Permanence Permanence -> Synapse
  ;; produce synapse with decremented permanence (trimmed to min)
  (let ((dec-perm (fx- (syn-perm synapse) by)))
    (make-syn (syn-prex synapse) (if (fx<? dec-perm trim)  min-perm  dec-perm))))
                                                                                            ;
(define (adapt-synapses                  ;; SP InputVec (listof ColumnX) ->
          sp input-vector active-columns)
  ;; update permanences in segments of active columns (+ if synapse's input on, - if not)
  (let ((syn-perm-trim-threshold (sp-syn-perm-trim-threshold sp))
        (syn-perm-active-inc     (sp-syn-perm-active-inc sp))
        (syn-perm-inactive-dec   (sp-syn-perm-inactive-dec sp))
        (syn-perm-connected      (sp-syn-perm-connected sp)))
    (for-each-segment sp (lambda (segment cx) ;; count connected & raise to threshold if needed
        (let loop ((i (fx1- (synapses-length segment))) (num-connected 0))
          (cond
            [ (fxnegative? i)
                (when (fx<? num-connected (sp-stimulus-threshold sp))
                  (raise-permanence-to-threshold sp segment)) ]
            [ else
              (let* (
                  [synapse (synapses-ref segment i)]
                  [synapse (if (bitwise-bit-set? input-vector (syn-prex synapse))
                               (increase-perm synapse syn-perm-active-inc)
                               (decrease-perm synapse syn-perm-inactive-dec syn-perm-trim-threshold))])
                (synapses-set! segment i synapse)
                (loop (fx1- i) (if (fx>=? (syn-perm synapse) syn-perm-connected)
                                   (fx1+ num-connected)
                                   num-connected)))])))
      active-columns)))
                                                                                            ;
(define (bump-up-weak-columns sp)        ;; SP ->
  ;; increase permanences of columns with odc below minodc
  (let ((weak-columns                    ;; (listof ColumnX) 
          (do ([cx (fx1- (sp-num-columns sp)) (fx1- cx)]
               [cs '() (if (fl<? (vector-ref (sp-overlap-duty-cycles sp) cx) 
                                 (vector-ref (sp-min-overlap-duty-cycles sp) cx)) 
                           (cons cx cs) 
                           cs)])
            ((fxnegative? cx ) cs)) ))
    (for-each-segment sp (lambda (segment cx)
        (do ([i (fx1- (synapses-length segment)) (fx1- i)]) ((fxnegative? i))
          (synapses-set! segment i (increase-perm (synapses-ref segment i) 
                                                (sp-syn-perm-below-stimulus-inc sp)))))
      weak-columns)))
                                                                                            ;
(define (raise-permanence-to-threshold   ;; SP Segment ->
          sp segment)
  ;; increase positive permanences of column until enough are connected
  (let* ( [syn-perm-connected (sp-syn-perm-connected sp)]
          [num-connected (synapses-count (lambda (synapse)
                             (fx>=? (syn-perm synapse) syn-perm-connected))
                           segment)])
    (when (fx<? num-connected (sp-stimulus-threshold sp))
      (let ((syn-perm-below-stimulus-inc (sp-syn-perm-below-stimulus-inc sp)))
        (do ((i (fx1- (synapses-length segment)) (fx1- i))) ((fxnegative? i))
          (when (positive? (synapses-ref segment i))
            (synapses-set! segment i 
                (increase-perm (synapses-ref segment i) syn-perm-below-stimulus-inc))))
        (raise-permanence-to-threshold sp segment)))))
                                                                                            ;
;; -- Initialize permanences --
  ;; see https://github.com/numenta/nupic/issues/3233
                                                                                            ;
#;(define (init-perm-connected sp)         ;; SP -> Permanence
  ;; produce random permanence above connected threshold
  (clip-max (fx+ (sp-syn-perm-connected sp)
     (fx3* (fx- max-perm (sp-syn-perm-connected sp)) (random max-perm)))))
                                                                                            ;
(define (init-perm-connected sp)         ;; SP -> Permanence
  ;; produce random permanence in small range above connected threshold
  (clip-max (fx+ (sp-syn-perm-connected sp)
                 (random (fx* 5 (sp-syn-perm-inactive-dec sp))))))
                                                                                            ;
#;(define (init-perm-non-connected sp)     ;; SP -> Permanence
  ;; produce random permanence below connected threshold
  (clip-min (fx3* (sp-syn-perm-connected sp) (random max-perm))))
                                                                                            ;
(define (init-perm-non-connected sp)     ;; SP -> Permanence
  ;; produce random permanence in small range below connected threshold
  (clip-min (fx- (sp-syn-perm-connected sp) 1
                 (random (fx* 5 (sp-syn-perm-active-inc sp))))))
                                                                                            ;
(define (init-permanence sp potential)   ;; SP (vectorof InputX) -> Segment
  ;; produce segment with initial permanences
  (build-synapses (vector-length potential) (lambda (i)
      (make-syn (vector-ref potential i)
        (if (fl<? (random 1.0) init-connected-pct)
            (init-perm-connected sp)
            (let ((ipnc (init-perm-non-connected sp)))
              (if (fx<? ipnc (sp-syn-perm-trim-threshold sp))  min-perm  ipnc)))))))
                                                                                            ;
(define (map-column sp cx)               ;; SP ColumnX -> InputX
  ;; produce index of input that will be centre of column's receptive field
  (let ((inp-dims (sp-input-dimensions sp))
        (col-dims (sp-column-dimensions sp)))
    (if (fx=? 1 (length inp-dims))
      (int<- (- (* (+ cx 0.5) (/ (sp-num-inputs sp) (sp-num-columns sp))) 0.5))
      (let* ( (column-coords (unravel2 cx col-dims))
              (ratios        (map / column-coords col-dims))
              (input-coords  (map exact (map floor (map + (map * inp-dims ratios) 
                                                          (map * (make-list (length inp-dims) 0.5)
                                                                 (map / inp-dims col-dims)))))))
        (ravel2 input-coords inp-dims)))))
                                                                                            ;
(define (map-potential sp cx)            ;; SP ColumnX -> (vectorof InputX)
  ;; produce vector of potential input indices for column
  (let* ((center-input (map-column sp cx))
         (column-inputs 
          (neighborhood center-input (sp-potential-radius sp) (sp-input-dimensions sp) (sp-wrap-around sp)))
         (num-potential (int<- (* (vector-length column-inputs) (sp-potential-pct sp)))))
    (vector-sample column-inputs num-potential)))
                                                                                            ;
;; -- Boost factors --
                                                                                            ;
(define (update-boost-factors sp)        ;; SP ->
  (sp-boost-factors-set! sp (if (sp-global-inhibition sp)
                                (update-boost-factors-global sp)
                                (update-boost-factors-local  sp))))
                                                                                            ;
(define (boost-func sp dc tgt-density)   ;; SP DutyCycle Number -> BoostFactor
  ;; produce boost factor as exponential of duty cycle
  (if (flzero? (sp-boost-strength sp)) 1.0
    (exp (fl- (fl* (sp-boost-strength sp) (fl- dc tgt-density))))))
                                                                                            ;
(define (update-boost-factors-global sp) ;; SP -> (ColVecOf BoostFactor)
  ;; produce boost factors for all columns from global target density
  (let ((target-density
          (if (flpositive? (sp-local-area-density sp)) 
              (sp-local-area-density sp)
              (let ((inhibition-area (fxmin (sp-num-columns sp)
                                            (expt (fx1+ (fx* 2 (sp-inhibition-radius sp)))
                                                  (length (sp-column-dimensions sp))))))
                (min 0.5 (/ (sp-num-active-columns-per-inh-area sp) (inexact inhibition-area)))))))
    (vector-map (lambda (adc)            ;; ActiveDutyCycle -> BoostFactor
        (boost-func sp adc target-density))
      (sp-active-duty-cycles sp))))
                                                                                            ;
(define (update-boost-factors-local sp)  ;; SP -> (ColVecOf BoostFactor)
  ;; produce boost factors for all columns from local target density
  (vector-map-x
    (lambda (adc cx)                     ;; DutyCycle ColumnX -> BoostFactor
      (let* ( (mask-neighbors 
                (neighborhood cx (sp-inhibition-radius sp) (sp-column-dimensions sp) (sp-wrap-around sp)))
              (target-density (vector-average (vector-refs (sp-active-duty-cycles sp) mask-neighbors))))
        (boost-func sp adc target-density)))
    (sp-active-duty-cycles sp)))
                                                                                            ;
(define (calculate-overlap sp input-vec) ;; SP InputVec -> (ColVecOf OverlapX)
  ;; produce overlaps: vector of counts of intersection of input and synapses
  (let ((connected-synapses (sp-connected-synapses sp)))
    (build-vector (sp-num-columns sp) (lambda (cx)
        (fixnum->flonum (bitwise-bit-count 
          (bitwise-and input-vec (vector-ref connected-synapses cx))))))))
                                                                                            ;
;; -- Inhibition --
                                                                                            ;
(define (inhibit-columns sp overlaps)    ;; SP (ColVecOf OverlapX) -> (listof ColumnX)
  ;; produce list of indices of active columns
  (let ((density
         (if (positive? (sp-local-area-density sp))
             (sp-local-area-density sp)
             (let ((inhibition-area (expt (fx1+ (fx* 2 (sp-inhibition-radius sp)))
                                          (length (sp-column-dimensions sp)))))
               (min 0.5 (/ (sp-num-active-columns-per-inh-area sp) 
                           (fxmin (sp-num-columns sp) inhibition-area)))))))
    (if (or (sp-global-inhibition sp)
            (fx>? (sp-inhibition-radius sp) (apply fxmax (sp-column-dimensions sp))))
        (inhibit-columns-global sp overlaps density)
        (inhibit-columns-local  sp overlaps density))))
                                                                                            ;
(define (inhibit-columns-global          ;; SP (ColVecOf OverlapX) Num -> (listof ColumnX)
          sp overlaps density)
  ;; produce indices of most active columns from given overlap counts and required density
  ;; when multiple columns have same overlap count, make random selection [cf .py "_tieBreaker"]
  (define  max-colx          (expt 2 (fxdiv (fixnum-width) 2)))
  (define (make-ovcx ov cx)  (fx+ (fx* (exact (flround (fl* ov 100000.0))) max-colx) cx))
  (define (ov ovcx)          (fxdiv ovcx max-colx))
  (define (cx ovcx)          (fxmod ovcx max-colx))
  (let ((ovcxs      (vector-map-x make-ovcx overlaps))
        (num-active (int<- (* density (sp-num-columns sp)))))
    (vector-sort! fx<? ovcxs)
    (let next ([swi (list)] [batch-top (fx1- (sp-num-columns sp))])   ;; swi = sorted-winner-indices
      (let ((batch-ov (ov (vector-ref ovcxs batch-top))))     ;; work down overlaps in batches with = count
        (let-values ([(this-batch next-top)
            (let batch ([swi-batch (list)] [bx batch-top])            ;; (listof ColumnX) ColumnX
              ;; produce batch and index to start next batch
              (let ((bx-ov (ov (vector-ref ovcxs bx))))
                (cond [(fxnegative? bx) (values swi-batch -1) ]
                      [(fx=? bx-ov batch-ov)
                       (batch (cons (cx (vector-ref ovcxs bx)) swi-batch) (fx1- bx))]
                      [else (values swi-batch bx)])))])
            (define (select)                                         ;; make random seln from last batch
              (let ((needed (fx- num-active (length swi))))
                (if (fxpositive? needed)
                  (vector->list (vector-sample (list->vector this-batch) needed))
                  (list) )))
          (cond [(fxnegative? next-top) (append swi (select))]       ;; no more overlaps
                [(fx<? (fx+ (length swi) (length this-batch)) num-active)
                 (next (append swi this-batch) next-top)]            ;; append & go for more
                [else (append swi (select))]))))))                   ;; enough active cols
                                                                                            ;
(define (tied-overlaps overlaps this-ov) ;; (vectorof OverlapX) (OverlapX) -> (vectorof ColumnX)
  ;; produce indices in neighborhood for overlaps with same count as this-ov
  (vector-filter id
    (vector-map-x (lambda (ov cx)
        (if (fl=? ov this-ov) cx #f))
      overlaps)))
                                                                                            ;
(define (inhibit-columns-local           ;; SP (ColVecOf OverlapX) Num -> (listof ColumnX)
          sp overlaps density)
  ;; produce indices of most active columns within each column's locality
  (let ((active (make-vector (sp-num-columns sp) #f)))
    (for-each-column sp
      (lambda (cx)                 ;; ColumnX -> [set active if in top [req density] cols in its neighborhood]
        (let ((this-ov (vector-ref overlaps cx)))
          (when (> this-ov (sp-stimulus-threshold sp))
            (let* ((neighbors      
                    (neighborhood cx (sp-inhibition-radius sp) (sp-column-dimensions sp) (sp-wrap-around sp)))
                   (nhood-overlaps (vector-refs overlaps neighbors))
                   (num-bigger     (vector-count (lambda (ov) (> ov this-ov)) nhood-overlaps))
                   (tied-neighbors (tied-overlaps nhood-overlaps this-ov))
                   (num-ties-lost  (vector-count (lambda (cx) (vector-ref active cx)) tied-neighbors)))
              (when (< (+ num-bigger num-ties-lost) (int<- (* density (vector-length neighbors))))
                (vector-set! active cx #t)))))))
    (vector-indices active)))
                                                                                            ;
;; -- Topology --
                                                                                            ;
(define (ravel2 coords dims)             ;; (listof Nat) (listof Nat) -> Nat
  ;; produce index for {row-coord, col-coord} in matrix of {nrows, ncols}
  (fx+ (fx* (car coords) (cadr dims)) (cadr coords)))
                                                                                            ;
(define (unravel2 index dims)            ;; Nat (listof Nat) -> (listof Nat)
  ;; produce {row-coord, col-coord} from index in matrix of {nrows,ncols}
  (let-values (((rows cols) (fxdiv-and-mod index (cadr dims))))
    (list rows cols)))
                                                                                            ;
(define (neighborhood center radius      ;; Nat Nat (listof Nat) Boolean -> (Vectorof Nat)
           dims wrap-around)
  ;; produce indices of elements within radius of center; 1 or 2 dimensional
  (let ((finish (lambda (center dim) (min (+ center radius 1) dim))))
    (if (null? (cdr dims))
      (if wrap-around                    ;; wrap-around only for 1 dimensional
        (let ([start (fxmod (fx- center radius) (car dims))]
              [n     (fx1+ (fx* 2 radius))])
          (build-vector n (lambda (i) (fxmod (fx+ start i) (car dims)))))
        (let* ( (start (fxmax 0 (fx- center radius)))
                (n (fx- (finish center (car dims)) start)))
          (build-vector n (lambda (i) (fx+ start i)))))
      (let* ( (center-coord (unravel2 center dims))                ;; 2 dimensional
              (start-row (fxmax 0 (- (car  center-coord) radius)))
              (start-col (fxmax 0 (- (cadr center-coord) radius)))
              (n-rows (fx- (finish (car  center-coord) (car  dims)) start-row))
              (n-cols (fx- (finish (cadr center-coord) (cadr dims)) start-col))
              (size (fx* n-rows n-cols))
              (v (make-vector size))
              (cx (ravel2 (list start-row start-col) dims))
              (skip (fx1+ (fx- (cadr dims) n-cols))))
        (do ((i 0 (fx1+ i))) ((fx=? i size) v)
          (vector-set! v i cx)
          (set! cx (if (fxzero? (fxmod (fx1+ i) n-cols))
                     (fx+ cx skip)
                     (fx1+ cx))))))))
                                                                                            ;
;; -- Connections --
                                                                                            ;
(define (connected-inputs sp cx)         ;; SP ColumnX -> InputVec
  ;; produce bitwise vector of whether permanence of synapses of column > connected threshold
  (let* ( (segment            (vector-ref (sp-potential-pools sp) cx))
          (syn-perm-connected (sp-syn-perm-connected sp)))
    (do ( [sx (fx1- (synapses-length segment)) (fx1- sx)]
          [connected 0
            (let ((synapse (synapses-ref segment sx)))
              (if (fx>=? (syn-perm synapse) syn-perm-connected)
                  (bitwise-copy-bit connected (syn-prex synapse) 1)
                  connected))])
        ((fxnegative? sx) connected))))
                                                                                            ;
(define (for-each-segment sp proc cols)  ;; SP (Segment ColumnX -> ) {ColumnX} ->
  ;; apply proc to segment of each col, and update connected caches for col
  (for-each (lambda (cx)
      (proc (vector-ref (sp-potential-pools sp) cx) cx)
      (let ((connected (connected-inputs sp cx)))
        (vector-set! (sp-connected-synapses sp) cx connected)
        (vector-set! (sp-connected-counts sp) cx (bitwise-bit-count connected))))
    cols))
                                                                                            ;
(define (for-each-column sp proc)  ;; SP (ColumnX -> ) ->
  ;; apply proc to each column index
  (do ([cx (fx1- (sp-num-columns sp)) (fx1- cx)]) ((fxnegative? cx ))
    (proc cx)))
                                                                                            ;
;; -- Smoke tests --
                                                                                            ;
(define (tests)
  (let ((SP42 (make-sp `(
                  [input-dimensions  . (4)]
                  [column-dimensions . (2)]
                  [global-inhibition . #t])))
        (COL (lambda ps 
              (build-synapses (length ps) (lambda (sx) (make-syn sx (list-ref ps sx)))))))

    [sp-potential-pools-set! SP42 (vector (COL 101) (COL 0 101 101))]
    [sp-connected-synapses-set! SP42
      (vector (connected-inputs SP42 0) (connected-inputs SP42 1))]

    [expect ( [connected-inputs SP42 0] #b0001 )  
            ( [connected-inputs SP42 1] #b0110 )]
    [expect ( [calculate-overlap SP42 #b1111] (vector 1.0 2.0) )
            ( [calculate-overlap SP42 #b1110] (vector 0.0 2.0) )
            ( [calculate-overlap SP42 #b0101] (vector 1.0 1.0) )]

    [sp-overlap-duty-cycles-set! SP42 '#(0 0)]
    [sp-active-duty-cycles-set!  SP42 '#(0 0)]

    [expect ((begin [update-duty-cycles SP42 '#(0.0 0.0) '()] 
              [sp-overlap-duty-cycles SP42]) '#(0.0 0.0) )]
    [expect ( [avg-columns-per-input SP42] 1/2 )]
    [expect ( [avg-connected-span-1d SP42 1] 2 )]

    (let ((SP2221 (make-sp `([input-dimensions . (2 2)] [column-dimensions . (2 1)]))))
      [sp-potential-pools-set! SP2221 (vector (COL 0 101 101 0) (COL 0))]
      [sp-connected-synapses-set! SP2221
        (vector (connected-inputs SP2221 0) (connected-inputs SP42 1))]
      [expect ( [avg-connected-span-2d SP2221 0] 1.5 )])
    )
    
  [expect ( [ravel2 '(1 1) '(99 2)] 3)]
  [expect ( [unravel2 3 '(99 2)] '(1 1) )]                                                                                        ;
  [expect ( [neighborhood  1 2 '(100) #f] '#(0 1 2 3)  )
          ( [neighborhood 99 2 '(100) #f] '#(97 98 99) )
          ( [neighborhood  0 2 '(5 5) #f] '#(0 1 2 5 6 7 10 11 12) )
          ( [neighborhood 12 2 '(5 5) #f] [build-vector 25 id] )
          ( [neighborhood 24 2 '(5 5) #f] '#(12 13 14 17 18 19 22 23 24) )
          ( [neighborhood 50 2 '(1 100) #f] '#(48 49 50 51 52) )]
  [expect ( [neighborhood  1 2 '(100) #t] '#(99 0 1 2 3)  )
          ( [neighborhood  3 2 '(100) #t] '#(1 2 3 4 5) )
          ( [neighborhood 99 2 '(100) #t] '#(97 98 99 0 1) )]

  (let ((vec5 (lambda (x) (make-vector 5 x))))
    [expect ( [update-duty-cycles-helper (vec5 1000) (vec5    0) 1000] (vec5  999) )
            ( [update-duty-cycles-helper (vec5 1000) (vec5 1000) 1000] (vec5 1000) )
            ( [update-duty-cycles-helper (vec5 1000) '#(2000 4000 5000 6000 7000) 1000]
                                                     '#(1001 1003 1004 1005 1006)  )])

  )
  
(tests)          
                                                                                            ;
)
