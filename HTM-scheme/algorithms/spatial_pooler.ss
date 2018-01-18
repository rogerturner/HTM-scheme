;; ========= HTM-scheme Spatial Pooler Copyright 2017 Roger Turner. =========
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
  ;; Plain Scheme using non-sparse vectors but packing values into Fixnums.
  ;; Max 2 dimensions, no wraparound, no parameter consistency checking.
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.

#!chezscheme
                                                                                            ;
(optimize-level 3)
                                                                                            ;
(library (spatial_pooler)
                                                                                            ;
(export
  perm<-
  overlap-count
  calculate-overlap
  make-sp*
  compute)
                                                                                            ;
(import 
  (rnrs)                   ;; use (except (chezscheme) add1 make-list random) for load-program
  (htm_prelude))
    
;; -- Spatial Pooler Types --
                                                                                            ;
;; X, Y, Z      = type parameters (arbitrary types in function type specification etc)
;; X Y -> Z     = function with argument types X, Y and result type Z
;; Nat          = natural number (including zero) (scheme Fixnum or exact Number)
;; Fix10k       = Fixnum interpreted as number with 4 decimal places, ie 10000 = 1.0
;; InputVec     = Number (Bignum), interpreted bitwise as (vectorof Boolean)
;; InputX       = Nat [0 - MAXINP], index of bit in InputVec
;; OverlapX     = Fixnum, interpreted as OverlapCount * MAXINP + InputX
;; Permanence   = Fix10k [0 - MAXPERM], interpreted as 0.0 - 1.0
;; Synapse      = Fixnum, interpreted as InputX * MAXPERM + Permanence 
;; Segment      = (vector of Synapse), the input indices and permanences of a column
;; DutyCycle    = Fix10k
;; BoostFactor  = Fix10k
;; ColumnX      = Nat [0 - MAXCOL], column index of Segment/DutyCycle/...
;; (ColVecOf X) = Vector of Xs indexed by ColumnX
;; Columns      = (ColVecOf Segment), (ColVecOf DutyCycle) ... the columns, part of SP
;; SP           = Record, sp parameters accessed as (sp-param-name sp)

;; -- Permanences and synapses --
                                                                                            ;
(define %max-perm% 9999)                 ;; Fixnum
(define %min-perm% 0)                    ;; Fixnum
                                                                                            ;
(define (perm<- x)                       ;; Number[0.0-.9999] -> Permanence
  (min %max-perm% (max %min-perm% (int<- (* x x10k)))))
                                                                                            ;
(define (make-synapse inpx perm)         ;; InputX Permanence -> Synapse
  (fx+ (fx* inpx x10k) perm))
                                                                                            ;
(define (inpx synapse)                   ;; Synapse -> InputX
  (fxdiv synapse x10k))
                                                                                            ;
(define (perm synapse)                   ;; Synapse -> Permanence
  (fxmod synapse x10k))
                                                                                            ;
;; -- Overlaps --
                                                                                            ;
(define %max-columns%                    ;; Fixnum
  (expt 2 (div (fixnum-width) 2)))
                                                                                            ;
(define (make-overlap count colx)        ;; Nat ColumnX -> OverlapX
  (fx+ (fx* count %max-columns%) colx))
                                                                                            ;
(define (count overlap)                  ;; OverlapX -> Nat
  (fxdiv overlap %max-columns%))
                                                                                            ;
(define (overlap-count overlap)          ;; OverlapX -> Nat
  (count overlap))
                                                                                            ;
(define (colx overlap)                   ;; OverlapX -> ColumnX
  (fxmod overlap %max-columns%))
                                                                                            ;
(define (list->overlaps l)               ;; (listof Count) -> (vectorof OverlapX)
  ;; (convenience for tests)
  (build-vector (length l) (lambda (i) (make-overlap (list-ref l i) i))))
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
(define (neighborhood center radius dims);; Nat Nat (listof Nat) -> (listof Nat)
  ;; produce indices of elements within radius of center; 1 or 2 dimensional
  (let ((finish (lambda (center dim) (min (+ center radius 1) dim))))
    (if (null? (cdr dims))
      (let* ( (start (fxmax 0 (fx- center radius)))                ;; 1 dimensional
              (n (fx- (finish center (car dims)) start)))
        (build-vector n (lambda (i) (fx+ start i))))
      (let* ( (center-coord (unravel2 center dims))                ;; 2 dimensional
              (start-row (fxmax 0 (- (car  center-coord) radius)))
              (start-col (fxmax 0 (- (cadr center-coord) radius)))
              (n-rows (fx- (finish (car  center-coord) (car  dims)) start-row))
              (n-cols (fx- (finish (cadr center-coord) (cadr dims)) start-col))
              (size (fx* n-rows n-cols))
              (v (make-vector size))
              (cx (ravel2 (list start-row start-col) dims))
              (skip (add1 (fx- (cadr dims) n-cols))))
        (do ((i 0 (add1 i))) ((fx=? i size) v)
          (vector-set! v i cx)
          (set! cx (if (fxzero? (fxmod (add1 i) n-cols))
                     (fx+ cx skip)
                     (add1 cx))))))))
                                                                                            ;
(define (for-each-column-index sp proc)  ;; SP (Nat -> ) ->
  ;; apply proc to each column index
  (do ((cx 0 (add1 cx))) ((fx=? cx (sp-num-columns sp)))
    (proc cx)))
                                                                                            ;
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
    (mutable boost-strength)                  ;; Number
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
    (mutable connected-synapses)))            ;; (ColVecOf InputVec)

  ;; cf https://github.com/numenta/nupic/pull/3143 and https://github.com/numenta/nupic/issues/3233
(define sp-defaults `(                   ;; (listof [Keyword . X])
    (potential-radius                . 16)
    (potential-pct                   . 0.5)
    (global-inhibition               . #f)
    (num-active-columns-per-inh-area . 10)
    (local-area-density              . -1.0)
    (stimulus-threshold              . 0)
    (inhibition-radius               . 0)
    (duty-cycle-period               . 1000)
    (boost-strength                  . 0.0)
    (iteration-num                   . 0)
    (iteration-learn-num             . 0)
    (update-period                   . 50)
    (syn-perm-trim-threshold         . ,(perm<-  0.025))
    (syn-perm-active-inc             . ,(perm<-  0.050))
    (syn-perm-inactive-dec           . ,(perm<-  0.008))
    (syn-perm-below-stimulus-inc     . ,(perm<-  0.010))
    (syn-perm-connected              . ,(perm<-  0.100))
    (min-pct-overlap-duty-cycles     . ,(fx10k<- 0.001))
    (boost-factors                   . ,'#())
    (overlap-duty-cycles             . ,'#())
    (active-duty-cycles              . ,'#())
    (min-overlap-duty-cycles         . ,'#())
    (potential-pools                 . ,'#())
    (connected-synapses              . ,'#())))
                                                                          ;
(define (sp-num-inputs  sp)              ;; SP -> Nat
  (apply * (sp-input-dimensions  sp)))
(define (sp-num-columns sp)              ;; SP -> Nat
  (apply * (sp-column-dimensions sp)))
(define init-connected-pct (fx10k<- 0.5));; Fix10k
                                                                                            ;
;; -- Connections --
                                                                                            ;
(define (connected-inputs sp cx)         ;; SP ColumnX -> InputVec
  ;; produce bitwise vector of whether permanence of synapses of column > connected threshold
  (let* ( (segment (vector-ref (sp-potential-pools sp) cx))
          (segment-length (vector-length segment))
          (syn-perm-connected (sp-syn-perm-connected sp)))
    (let next ((synx 0) (connected 0))
      (cond [(fx=? synx segment-length) connected]
            [else (next (add1 synx)
                    (let ((synapse (vector-ref segment synx)))
                      (if (fx>? (perm synapse) syn-perm-connected)
                          (bitwise-copy-bit connected (inpx synapse) 1)
                          connected)))]))))
                                                                                            ;
(define (update-connected-synapses sp cx);; SP ColumnX ->
  ;; update entry in the connected synapses cache
  (vector-set! (sp-connected-synapses sp) cx (connected-inputs sp cx)))
                                                                                            ;
(define (connected-synapses sp)          ;; SP -> (vectorof InputVec)
  ;; create connected synapses cache
  (build-vector (sp-num-columns sp) (lambda (cx) (connected-inputs sp cx))))
                                                                                            ;
(define (default-sp id cd . args)        ;; (listof Nat) (listof Nat) (listof KWarg) -> SP
  ;; produce spatial pooler with default/dummy parameters for use in tests
  (let ((sp (apply make-sp (append (list id cd) (key-word-args args sp-defaults)))))
    (when (zero? (vector-length (sp-potential-pools sp)))
      (sp-potential-pools-set! sp (make-vector (sp-num-columns sp) '#())))
    (sp-connected-synapses-set! sp (connected-synapses sp))
    sp))
                                                                                            ;
(define (boosted-overlaps sp overlaps)   ;; SP (ColVecOf OverlapX) -> (ColVecOf OverlapX)
  ;; produce overlaps adjusted by boost factors
  (vector-map
    (lambda (ov bf)
      (make-overlap (int<- (fx10k* bf (count ov))) (colx ov)))
    overlaps (sp-boost-factors sp)))
                                                                                            ;
;; -- Duty cycles --
                                                                                            ;
(define (update-min-duty-cycles sp)      ;; SP ->
  ;; recompute minimum overlap duty cycles as percentage of maximum
  (if (or (sp-global-inhibition sp) (> (sp-inhibition-radius sp) (sp-num-inputs sp)))
      (vector-fill! (sp-min-overlap-duty-cycles sp) (min-duty-cycles-global sp))  ;; global
      (for-each-column-index sp                                                   ;; local
        (lambda (cx)
          (vector-set! (sp-min-overlap-duty-cycles sp) cx (min-duty-cycles-local sp cx))))))
                                                                                            ;
(define (min-duty-cycles-global sp)      ;; SP -> DutyCycle
  ;; produce minodc as percent of max odc of all columns
  (fx10k* (sp-min-pct-overlap-duty-cycles sp) (fxvector-max (sp-overlap-duty-cycles sp))))
                                                                                            ;
(define (min-duty-cycles-local sp cx)    ;; SP ColumnX -> DutyCycle
  ;; produce minodc of column as percent of max odc of columns in its neighborhood
  (let* ( (locality         (neighborhood cx (sp-inhibition-radius sp) (sp-column-dimensions sp)))
          (max-overlap-duty (fxvector-max (vector-refs (sp-overlap-duty-cycles sp) locality))))
    (fx10k* (sp-min-pct-overlap-duty-cycles sp) max-overlap-duty)))
                                                                                            ;
(define (update-duty-cycles-helper       ;; (ColVecOf DutyCycle) (ColVecOf Number) Nat -> (ColVecOf DutyCycle)
          duty-cycles new-input period)  
  ;; produce updated duty cycle estimate
  (vector-map (lambda (dc ni)
                (int<- (/ (+ (* dc (- period 1)) ni) period)))
              duty-cycles new-input))
                                                                                            ;
(define (update-duty-cycles              ;; SP (ColVecOf OverlapX) (listof ColumnX) ->
          sp overlaps active-columns)
  ;; recompute duty cycles from this input's overlaps and resulting active columns
  (let* ((overlap-array (vector-map (lambda (ov) (if (positive? (count ov)) x10k 0))
                                    overlaps))
         (active-array  (make-vector (sp-num-columns sp) 0)))
    (for-each 
      (lambda (cx) (vector-set! active-array cx x10k)) 
      active-columns)
    (sp-overlap-duty-cycles-set! sp 
      (update-duty-cycles-helper (sp-overlap-duty-cycles sp) overlap-array (sp-duty-cycle-period sp)))
    (sp-active-duty-cycles-set!  sp 
      (update-duty-cycles-helper (sp-active-duty-cycles  sp) active-array  (sp-duty-cycle-period sp)))))
                                                                                            ;
;; -- Inhibition radius --
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
      (cond [(= rx nrows) (/ (+ (apply max l) (- (length l) 1)) 2)]
            [else
              (let* ((rowstart (* rx ncols))
                     (row (bitwise-bit-field connected rowstart (+ rowstart ncols))))
                (next-row (add1 rx) 
                          (if (zero? row) l
                            (cons (bitwise-span row) l))))]))))
                                                                                            ;
(define (inhibition-radius sp)           ;; SP -> Nat
  ;; produce new inhibition radius
  (if (sp-global-inhibition sp)
      (apply max (sp-column-dimensions sp))         ;; global
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
;; -- Adapt synapses --
                                                                                            ;
(define (for-each-col sp proc cols)      ;; SP (Segment -> ) (listof ColumnX) ->
  ;; apply proc to segment of each column in list, and update connected cache for col
  (for-each
    (lambda (cx)
      (proc (vector-ref (sp-potential-pools sp) cx))
      (update-connected-synapses sp cx))
    cols))
                                                                                            ;
(define (increase-perm synapse by)       ;; Synapse Permanence -> Synapse
  ;; produce synapse with incremented permanence (clipped to max)
  (make-synapse (inpx synapse) (fxmin (fx+ (perm synapse) by) %max-perm%)))
                                                                                            ;
(define (decrease-perm synapse by trim)  ;; Synapse Permanence Permanence -> Synapse
  ;; produce synapse with decremented permanence (trimmed to min)
  (let ((dec-perm (fx- (perm synapse) by)))
    (make-synapse (inpx synapse) (if (fx<? dec-perm trim)
                                     %min-perm%
                                     dec-perm))))
                                                                                            ;
(define (raise-permanence-to-threshold   ;; SP Segment ->
          sp segment)
  ;; increase all permanences of column until enough are connected
  (let ((syn-perm-below-stimulus-inc (sp-syn-perm-below-stimulus-inc sp))
        (syn-perm-connected (sp-syn-perm-connected sp)))
    (do ((i (fx- (vector-length segment) 1) (fx- i 1))) ((fxnegative? i))
      (vector-set! segment i (increase-perm (vector-ref segment i) syn-perm-below-stimulus-inc)))
    (let ((num-connected (vector-count 
                           (lambda (synapse)
                             (fx>=? (perm synapse) syn-perm-connected))
                           segment)))
      (when (fx<? num-connected (sp-stimulus-threshold sp))
        (raise-permanence-to-threshold sp segment)))))
                                                                                            ;
(define (adapt-synapses                  ;; SP InputVec (listof ColumnX) ->
          sp input-vector active-columns)
  ;; update permanences in segments of active columns (+ if synapse's input on, - if not)
  (let ((syn-perm-trim-threshold (sp-syn-perm-trim-threshold sp))
        (syn-perm-active-inc     (sp-syn-perm-active-inc sp))
        (syn-perm-inactive-dec   (sp-syn-perm-inactive-dec sp))
        (syn-perm-connected      (sp-syn-perm-connected sp)))
    (for-each-col sp
      (lambda (segment)                  ;; count connected & raise to threshold if needed
        (let loop ((i (fx- (vector-length segment) 1)) (num-connected 0))
          (cond [(fxnegative? i) (when (fx<? num-connected (sp-stimulus-threshold sp))
                                       (raise-permanence-to-threshold sp segment))]
                [else
                  (let* ( (synapse (vector-ref segment i))
                          (synapse (if (bitwise-bit-set? input-vector (inpx synapse))
                                       (increase-perm synapse syn-perm-active-inc)
                                       (decrease-perm synapse syn-perm-inactive-dec syn-perm-trim-threshold))))
                    (vector-set! segment i synapse)
                    (loop (fx- i 1) (if (fx>=? (perm synapse) syn-perm-connected)
                                        (add1 num-connected)
                                        num-connected)))])))
      active-columns)))
                                                                                            ;
(define (bump-up-weak-columns sp)        ;; SP ->
  ;; increase permanences of columns with odc below minodc
  (let ((weak-columns                    ;; (listof ColumnX) 
         (do ((cx 0 (add1 cx))
              (l '() (if (< (vector-ref (sp-overlap-duty-cycles sp) cx) 
                            (vector-ref (sp-min-overlap-duty-cycles sp) cx)) 
                         (cons cx l) 
                         l)))
             ((= cx (sp-num-columns sp)) l))))
    (for-each-col sp
      (lambda (segment)
        (do ((i (fx- (vector-length segment) 1) (fx- i 1))) ((fxnegative? i))
          (vector-set! segment i (increase-perm (vector-ref segment i) 
                                                (sp-syn-perm-below-stimulus-inc sp)))))
      weak-columns)))
                                                                                            ;
;; -- Initialize permanences --
                                                                                            ;
(define (init-perm-connected sp)         ;; SP -> Permanence
  ;; produce random permanence in small range above connected threshold
  (+ (sp-syn-perm-connected sp) (random (* 10 (sp-syn-perm-inactive-dec sp)))))
                                                                                            ;
(define (init-perm-non-connected sp)     ;; SP -> Permanence
  ;; produce random permanence in small range below connected threshold
  (- (sp-syn-perm-connected sp) (random (* 10 (sp-syn-perm-active-inc sp)))))
                                                                                            ;
(define (init-permanence sp potential)   ;; SP (vectorof InputX) -> Segment
  ;; produce segment with initial permanences
  (vector-map
    (lambda (input-index)
      (make-synapse input-index
                    (if (< (random 10000) init-connected-pct)
                        (init-perm-connected sp)
                        (init-perm-non-connected sp))))
    potential))
                                                                                            ;
(define (map-column sp cx)               ;; SP ColumnX -> InputX
  ;; produce index of input that will be centre of column's receptive field
  (let ((inp-dims (sp-input-dimensions sp))
        (col-dims (sp-column-dimensions sp)))
    (if (= 1 (length inp-dims))
      (int<- (- (* (+ cx 0.5) (/ (sp-num-inputs sp) (sp-num-columns sp))) 0.5))
      (let* ( (column-coords (unravel2 cx col-dims))
              (ratios        (map / column-coords col-dims))
              (input-coords  (map exact (map floor (map + (map * inp-dims ratios) 
                                                          (map * (make-list (length inp-dims) 0.5)
                                                                 (map / inp-dims col-dims)))))))
        (ravel2 input-coords inp-dims)))))
                                                                                            ;
(define (map-potential-v sp cx)          ;; SP ColumnX -> (vectorof InputX)
  ;; produce vector of potential input indices for column
  (let* ((center-input (map-column sp cx))
         (column-inputs (neighborhood center-input (sp-potential-radius sp) (sp-input-dimensions sp)))
         (num-potential (int<- (* (vector-length column-inputs) (sp-potential-pct sp)))))
    (vector-sample column-inputs num-potential)))
                                                                                            ;
;; -- Boost factors --
                                                                                            ;
(define (update-boost-factors sp)        ;; SP ->
  (sp-boost-factors-set! sp (if (sp-global-inhibition sp)
                                (boost-factors-global sp)
                                (boost-factors-local  sp))))
                                                                                            ;
(define (boost-func sp dc tgt-density)   ;; SP DutyCycle Number -> BoostFactor
  ;; produce boost factor as exponential of duty cycle
  (if (zero? (sp-boost-strength sp)) x10k
    (fx10k<- (exp (- (* (sp-boost-strength sp) (- (/ dc x10k) tgt-density)))))))
                                                                                            ;
(define (boost-factors-global sp)        ;; SP -> (ColVecOf BoostFactor)
  ;; produce boost factors for all columns from global target density
  (let ((target-density
          (if (positive? (sp-local-area-density sp)) (sp-local-area-density sp)
              (let ((inhibition-area (min (sp-num-columns sp)
                                            (expt (add1 (* 2 (sp-inhibition-radius sp)))
                                                  (length (sp-column-dimensions sp))))))
                (min 0.5 (/ (sp-num-active-columns-per-inh-area sp) (inexact inhibition-area)))))))
    (vector-map
      (lambda (adc)                      ;; DutyCycle -> BoostFactor
        (boost-func sp adc target-density))
      (sp-active-duty-cycles sp))))
                                                                                            ;
(define (boost-factors-local sp)         ;; SP -> (ColVecOf BoostFactor)
  ;; produce boost factors for all columns from local target density
  (vector-map-indexed
    (lambda (adc cx)                     ;; DutyCycle ColumnX -> BoostFactor
      (let* ( (mask-neighbors (neighborhood cx (sp-inhibition-radius sp) (sp-column-dimensions sp)))
              (target-density (/ (vector-average (vector-refs (sp-active-duty-cycles sp) mask-neighbors)) 10000.)))
        (boost-func sp adc target-density)))
    (sp-active-duty-cycles sp)))
                                                                                            ;
(define (calculate-overlap sp input-vec) ;; SP InputVec -> (ColVecOf OverlapX)
  ;; produce overlaps: vector of counts of intersection of input and synapses
  (build-vector 
    (sp-num-columns sp)
    (lambda (cx)
      (let ((count (bitwise-bit-count 
              (bitwise-and input-vec (vector-ref (sp-connected-synapses sp) cx)))))
        (make-overlap count cx)))))
                                                                                            ;
;; -- Inhibition --
                                                                                            ;
(define (inhibit-columns-global          ;; SP (ColVecOf OverlapX) Num -> (listof ColumnX)
          sp overlaps density)
  ;; produce indices of most active columns from given overlap counts and required density
  ;; when multiple columns have same overlap count, make random selection [cf .py "_tieBreaker"]
  (vector-sort! [lambda (x y) (< (count x) (count y))] overlaps)
  (let ((num-active (int<- (* density (sp-num-columns sp)))))
    (let next ((swi '()) (this-top (- (sp-num-columns sp) 1)))         ;; swi = sorted-winner-indices
      (let ((this-batch-count (count (vector-ref overlaps this-top)))) ;; work down overlaps in batches with = count
        (let-values ([(this-batch next-top)
          (let batch ((swi-batch '()) (bx this-top))  ;; (listof ColumnX) ColumnX -> (listof ColumnX) ColumnX
            ;; produce batch and index to start next batch
            (let ((bx-count (count (vector-ref overlaps bx))))
              (cond [(or (zero? bx) (< bx-count (sp-stimulus-threshold sp)))
                     (values swi-batch -1)]
                    [(= bx-count this-batch-count)
                     (batch (cons (colx (vector-ref overlaps bx)) swi-batch) (- bx 1))]
                    [else (values swi-batch bx)])))])
          (let ((last-batch (lambda ()                                 ;; fn to make random seln from last batch
                  (let ((needed (- num-active (length swi))))
                    (if (positive? needed)
                      (vector->list (vector-sample (list->vector this-batch) needed))
                      '() )))))
            (cond [(negative? next-top) (append swi (last-batch))]     ;; no more overlaps
                  [(<= (+ (length swi) (length this-batch)) num-active)       
                   (next (append swi this-batch) next-top)]            ;; append & go for more
                  [else (append swi (last-batch))])))))))              ;; enough active cols
                                                                                            ;
(define (tied-overlaps overlaps this-ov) ;; (vectorof OverlapX) (OverlapX) -> (vectorof ColumnX)
  ;; produce indices in neighborhood for overlaps with same count as this-ov
  (vector-filter id
    (vector-map-indexed
      (lambda (ov cx)
        (if (= (count ov) (count this-ov)) cx #f))
      overlaps)))
                                                                                            ;
(define (inhibit-columns-local           ;; SP (ColVecOf OverlapX) Num -> (listof ColumnX)
          sp overlaps density)
  ;; produce indices of most active columns within each column's locality
  (let ((active (make-vector (sp-num-columns sp) #f)))
    (for-each-column-index sp
      (lambda (cx)                 ;; ColumnX -> [set active if in top [req density] cols in its neighborhood]
        (let ((this-ov (vector-ref overlaps cx)))
          (when (> (count this-ov) (sp-stimulus-threshold sp))
            (let* ((neighbors      (neighborhood cx (sp-inhibition-radius sp) (sp-column-dimensions sp)))
                   (nhood-overlaps (vector-refs overlaps neighbors))
                   (num-bigger     (vector-count (lambda (ov) (> ov this-ov)) nhood-overlaps))
                   (tied-neighbors (tied-overlaps nhood-overlaps this-ov))
                   (ties           (vector-refs nhood-overlaps tied-neighbors))
                   (num-ties-lost  (vector-count (lambda (ov) (vector-ref active (colx ov))) ties)))
              (when (< (+ num-bigger num-ties-lost) (int<- (* density (vector-length neighbors))))
                (vector-set! active cx #t)))))))
    (vector-indices active)))
                                                                                            ;
(define (inhibit-columns sp overlaps)    ;; SP (ColVecOf OverlapX) -> (listof ColumnX)
  ;; produce list of indices of active columns
  (let ((density
         (if (positive? (sp-local-area-density sp))
             (sp-local-area-density sp)
             (let ((inhibition-area (expt (add1 (* 2 (sp-inhibition-radius sp)))
                                          (length (sp-column-dimensions sp)))))
               (min 0.5 (/ (sp-num-active-columns-per-inh-area sp) 
                           (min (sp-num-columns sp) inhibition-area)))))))
    (if (or (sp-global-inhibition sp)
            (> (sp-inhibition-radius sp) (apply max (sp-column-dimensions sp))))
        (inhibit-columns-global sp overlaps density)
        (inhibit-columns-local  sp overlaps density))))
                                                                                            ;
;; -- Initialization --
                                                                                            ;
(define (make-sp* id cd . args)          ;; (listof Nat) (listof Nat) (listof KWarg) -> SP
  ;; produce spatial pooler with inp/col dims, default parameters overridden by args
  (let ((sp (apply default-sp (append (list id) (list cd) args))))
    (sp-syn-perm-trim-threshold-set!     sp (div (sp-syn-perm-active-inc sp) 2))
    (sp-syn-perm-below-stimulus-inc-set! sp (div (sp-syn-perm-connected sp) 10))
    (sp-potential-pools-set!             sp (build-vector (sp-num-columns sp)
                                              (lambda (cx)
                                                (init-permanence sp (map-potential-v sp cx)))))
    (sp-boost-factors-set!               sp (make-vector (sp-num-columns sp) x10k))
    (sp-overlap-duty-cycles-set!         sp (make-vector (sp-num-columns sp) 0))
    (sp-active-duty-cycles-set!          sp (make-vector (sp-num-columns sp) 0))
    (sp-min-overlap-duty-cycles-set!     sp (make-vector (sp-num-columns sp) 0))
    (sp-connected-synapses-set!          sp (connected-synapses sp))
    (sp-inhibition-radius-set!           sp (inhibition-radius sp))
    sp))
                                                                                            ;
(define (compute sp input-vector learn)  ;; SP InputVec Boolean -> (listof ColumnX)
  ;; produce active columns from input; optionally update sp if learning
  (sp-iteration-num-set! sp (add1 (sp-iteration-num sp)))
  (when learn (sp-iteration-learn-num-set! sp (add1 (sp-iteration-learn-num sp))))
  (let* ((overlaps (calculate-overlap sp input-vector))         ;; (ColVecOf OverlapX)
         (active-columns                                        ;; (listof ColumnX)
           (inhibit-columns sp (if learn
                                   (boosted-overlaps sp overlaps)
                                   overlaps))))
    (when learn
      (adapt-synapses sp input-vector active-columns)
      (update-duty-cycles sp overlaps active-columns)
      (bump-up-weak-columns sp)
      (update-boost-factors sp)
      (when (zero? (mod (sp-iteration-num sp) (sp-update-period sp)))
          (sp-inhibition-radius-set! sp (inhibition-radius sp))
          (update-min-duty-cycles sp)))
    active-columns))
                                                                                            ;
)