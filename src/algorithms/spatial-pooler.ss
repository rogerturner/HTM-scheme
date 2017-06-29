#!r6rs
;; ========= HTM-scheme Spatial Pooler Copyright 2017 Roger Turner. =========
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Based on spatial_pooler.py which is part of the Numenta Platform for  ;;
  ;; Intelligent Computing (NuPIC) Copyright (C) 2013-2016, Numenta, Inc.  ;;
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

  ;; Scheme top-level program with some tests and hello_sp and sp_tutorial examples.
  ;; Run in DrRacket or load in Chez Scheme repl then (hello-sp) or (sp-tutorial)

  ;; Translated from NuPIC spatial_pooler.py, see comments there for more info.
  ;; Plain Scheme using non-sparse vectors but packing values into Fixnums.
  ;; Max 2 dimensions, no wraparound, no parameter consistency checking.

  ;; See "(define (list-average..." for example of type definition comment and test.
  ;; Use a "Fold All" view (in eg Atom) for a source overview.
  
  (import (rnrs))  ;; use (except (chezscheme) add1 make-list random) for Chez load-program
  
  ;; "The last thing one discovers in composing a work is what to put first" (Pascal)

;; -- Spatial Pooler Types --
                                                                                            ;
;; X, Y, Z      = type parameters (arbitrary types in function type specification etc)
;; X Y -> Z     = function with argument types X, Y and result type Z
;; Nat          = natural number (including zero) (Scheme Fixnum or exact Integer)
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

;; ================== Prelude: convenience functions/style examples ==================

(define-syntax expect                    ;; ((X ... -> Y) X ...) Y -> [error]
  ;; check that function application(s) to arguments match expected values
  (lambda (x)                            
    (syntax-case x ()                                  ;; [expect ([fn args] expected ) ... ]
      [ (_ (expr expected) ...)                        ;; expr matches [fn args]
        #'(begin (let ((result expr))                  ;; eval expr just once
                   (unless (equal? result expected)    ;; no output if check passes
                           (for-each display `("**" expr #\newline 
                             "  expected: " expected #\newline 
                             "  returned: " ,result  #\newline)))) ...)])))
                                                                                            ;
(define (list-average l)                 ;; (listof Number) -> Number
  ;; produce mean of non-empty list <-      -------------------------
  (/ (apply + l) (length l)))        ;            ^
                                     ;            |
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; | ;;;;;;;;;; | ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; documentation for function: Description, Type, and Test (*)                ;;
      ;;                                                    |                       ;;
  [expect ( [list-average (list 1 2 27)]  10 )] ;; <---------                       ;;
      ;;    ----------------------------  --                                        ;;
      ;;    ^ exercise function           ^ expected result                         ;;
      ;; (*)  description and test are optional                                     ;;
      ;;      locate uses of fn by finding "(fn ", tests use "[fn ...] ")           ;;
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;                                                                                 
                                                                                            ;
(define (add1 n)                         ;; Fixnum -> Fixnum
  (fx+ 1 n))
                                                                                            ;
(define (build-list n proc)              ;; Nat (Nat -> X) -> (listof X)
  ;; produce list of n X's by applying proc to 0..n-1
  (let loop [ (i 0) (xs '()) ]
    (cond [(fx=? i n) (reverse xs)]
          [else (loop (add1 i) (cons (proc i) xs))])))
          
  [expect ( [build-list 3 -] '(0 -1 -2) )]
                                                                                            ;
(define (take n xs)                      ;; Nat (listof X) -> (listof X)
  ;; produce first n elements of xs
  (let loop [ (n n) (xs xs) (ys '()) ]
    (cond [(or (fxzero? n) (null? xs)) (reverse ys)]
          [else (loop (fx- n 1) (cdr xs) (cons (car xs) ys))])))
          
  [expect ( [take 2 '(a 2 c)] '(a 2) )]
                                                                                            ;
(define (make-list n x)                  ;; Nat X -> (listof X)
  ;; produce list of n x's
  (vector->list (make-vector n x)))

  [expect ( [make-list 0 9] '()      )
          ( [make-list 3 9] '(9 9 9) )]
                                                                                            ;
(define (list->bitwise ns)               ;; (listof Nat) -> Number
  ;; produce bitwise value by setting bits indexed by elements of ns
  (fold-left (lambda (acc n) (bitwise-copy-bit acc n 1)) 0 ns))
  
  [expect ( [list->bitwise '(0 2 3)] #b1101 )]
                                                                                            ;
(define (int<- x)                        ;; Number -> Integer
  (exact (round x)))
  
  [expect ( [int<- 2.500000000000001] 3 )
          ( [int<- 5/2              ] 2 )
          ( [int<- 7/2              ] 4 )]
                                                                                            ;
(define (id x)                           ;; X -> X
  ;; produce argument (convenient for tests)
  x)
                                                                                            ;
(define (build-vector n proc)            ;; Nat (Nat -> X) -> (vectorof X)
  ;; produce vector length n by applying proc to indices
  (let ((vec (make-vector n)))
    (do ((i 0 (add1 i))) ((fx=? i n) vec)
      (vector-set! vec i (proc i)))))
      
  [expect ( [build-vector 3 id] '#(0 1 2) )]
                                                                                            ;
(define (vector-filter pred vec)         ;; (X -> Bool) (vectorof X) -> (vectorof X)
  ;; produce vector of elements of vec for which (pred elt) is not false
  (list->vector (filter pred (vector->list vec))))

  [expect ( [vector-filter even? '#(1 2 3 4 5)] '#(2 4) )]
                                                                                            ;
(define (vector-map-indexed proc vec)    ;; (X Nat -> Y) (vectorof X) -> (vectorof Y)
  ;; produce vector by applying proc to each element of vec and its index
  (let ((result (make-vector (vector-length vec))))
    (do ((i 0 (add1 i))) ((fx=? i (vector-length vec)) result)
      (vector-set! result i (proc (vector-ref vec i) i)))))
    
  [expect ( [vector-map-indexed (lambda (x i) (+ x i)) '#(10 11 12)] '#(10 12 14) )]
                                                                                            ;
(define (vector-fold-left proc obj vec)  ;; (X Y -> X) X (vectorof Y) -> X
  ;; (proc ... (proc (proc obj vec-0) vec-1) ... vec-last)
  (let ((acc obj))
    (vector-for-each
      (lambda (x) (set! acc (proc acc x)))
      vec)
    acc))

  [expect ( [vector-fold-left (lambda (l x) (cons x l)) '() '#(0 1 2)] '(2 1 0) )]
                                                                                            ;
(define (vector-average vec)             ;; (vectorof Number) -> Number
  ;; (list-average (vector->list vec))
  (/ (vector-fold-left + 0 vec) (vector-length vec)))

  [expect ( [vector-average '#(1 2 5)] 8/3 )]
                                                                                            ;
(define (fxvector-max vec)               ;; (vectorof Fixnum) -> Fixnum
  ;; produce max element of vec
  (vector-fold-left fxmax (least-fixnum) vec))

  [expect ( [fxvector-max '#(0 -3 7 5)] 7 )]
                                                                                            ;
(define x10k 10000)                      ;; Fixnum                                                                                            
                                                                                            ;
(define (fx10k<- x)                      ;; [vectorof] Number -> [vectorof] Fix10k
  ;; produce Fix10k representation of number or vector of numbers
  (if (vector? x) (vector-map fx10k<- x)
                  (int<- (* x x10k))))
                  
  [expect ( [fx10k<- '#(1/3 1.0)] '#(3333 10000) )]
                                                                                            ;
(define (fx10k* x y)                     ;; Fix10k Fix10k -> Fix10k
  ;; produce product
  (div (+ (* x y) (div x10k 2)) x10k))
  
  [expect ( [fx10k* 5000 5000] 2500 )]
                                                                                            ;
(define (vector-take vec n)              ;; (vectorof X) Nat -> (vectorof X)
  ;; produce copy of first n elements of vec, or copy of vec if n > length
  (let* ( (size (fxmin n (vector-length vec)))
          (result (make-vector size)))
    (do ((i 0 (add1 i))) ((fx=? i size) result)
      (vector-set! result i (vector-ref vec i)))))

  [expect ( [vector-take '#(0 1 2 3) 0] '#()       )
          ( [vector-take '#(0 1 2 3) 2] '#(0 1)    )
          ( [vector-take '#(0 1 2 3) 5] '#(0 1 2 3))]
                                                                                            ;
(define (vector-count proc vec)          ;; (X -> Bool) (vectorof X) -> Nat
  ;; produce count of elements of vec for which (proc elt) is not false
  (vector-fold-left 
    (lambda (acc elt) (if (proc elt) (add1 acc) acc))
    0 vec))
    
  [expect ( [vector-count zero? '#(0 1 2 0 3)] 2 )
          ( [vector-count not   '#(#t 0 add1)] 0 )]
                                                                                            ;
(define (vector-indices vec)             ;; (vectorof X) -> (listof Nat)
  ;; produce list of indices for which vector element is not false
  (do [ (i 0 (add1 i)) 
        (l '() (if (vector-ref vec i) (cons i l) l))]
      ((fx=? i (vector-length vec)) l )))

  [expect ( [vector-indices '#(#f #f #f)] '()    )
          ( [vector-indices '#(#t #f  0)] '(2 0) )]
                                                                                            ;
(define (vector-refs vec refs)           ;; (vectorof X) (vectorof Nat) -> (vectorof X)
  ;; produce vector of selected elements of vec indexed by refs
  (let ((vrefs (make-vector (vector-length refs))))
    (do ((i 0 (add1 i))) ((fx=? i (vector-length refs)) vrefs)
      (vector-set! vrefs i (vector-ref vec (vector-ref refs i))))))

  [expect ( [vector-refs '#(0 11 22 33 44) '#(0 2 4)] '#(0 22 44) )]
                                                                                            ;
(define (vector->bitwise vec)            ;; (vectorof Nat) -> Number
  (list->bitwise (vector->list vec)))
                                                                                            ;
(define (bitwise-span n)                 ;; Number -> Nat
  ;; produce span of 1 bits in n
  (fx- (bitwise-length n) (bitwise-first-bit-set n)))
  
  [expect ( [bitwise-span #b01011101000] 7 )]
                                                                                            ;
(define (random n)                       ;; Nat -> Nat
  ;; produce random integer in range 0..n-1
  (let ((Lehmer-modulus   2147483647)
        (Lehmer-multiplier     48271))
    (set! *random-state* (mod (* *random-state* Lehmer-multiplier) Lehmer-modulus))
    (mod *random-state* n)))

  (define *random-state* 48271)
                                                                                            ;
(define (vector-sample source size)      ;; (vectorof X) Nat -> (vectorof X)
  ;; produce random selection of length size from source (using Durstenfeld shuffle)
  (let* ( (source-length (vector-length source))
          (shuffle (if (fx<? size source-length) size (fx- size 1))))
    (cond [(fxzero? size) '#() ]
          [(fx=? size 1)
           (vector (vector-ref source (random source-length))) ]
          [else
            (let ((vec (vector-take source source-length)))
              (do ((n 0 (add1 n))) ((fx=? n shuffle) (vector-take vec size))
                (let* ((r (fx+ n (random (fx- source-length n))))
                       (t (vector-ref vec r)))
                  (vector-set! vec r (vector-ref vec n))
                  (vector-set! vec n t))))])))

  [let* ( (source (build-vector 100 (lambda (n) n)))
          (select (vector->list [vector-sample source 50] )))
    (expect ( (length select) 50 )
            ( (list? (for-all memq select (make-list 50 (vector->list source)))) #t )  ;; all from source
            ( (map (lambda (x) (length (remq x select))) select) (make-list 50 49)) )] ;; all different
                                                                                            ;
(define (key-word-args args defaults)    ;; (listof KWarg) (listof KWarg) -> (listof X)
  ;; KWarg is [key . value]; produce list of default values overridden by arg values
  (for-each (lambda (arg-kv)
              (unless (assq (car arg-kv) defaults)
                (error #f "unknown keyword arg" arg-kv)))
            args)
  (map (lambda (default-kv)
         (let ((kv (assq (car default-kv) args)))  ;; if this default in args
           (if kv (cdr kv) (cdr default-kv))))     ;; then use given val
       defaults))

  [expect ( [(lambda (x . args)          ;; example fn with kws k1 & k2: default is (x 1 2)
              (append (list x) [key-word-args args '([k1 . 1] [k2 . 2])] ))  ;; end of fn defn
            99 '[k2 . 22] '[k1 . 11] #;'[k3 . 33] ] '(99 11 22) )]  ;; apply fn to 99 '[k2 ...

;; ================================= Spatial Pooler =================================

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
  
  [expect ( [ravel2 '(1 1) '(99 2)] 3)]
                                                                                            ;
(define (unravel2 index dims)            ;; Nat (listof Nat) -> (listof Nat)
  ;; produce {row-coord, col-coord} from index in matrix of {nrows,ncols}
  (let-values (((rows cols) (fxdiv-and-mod index (cadr dims))))
    (list rows cols)))
  
  [expect ( [unravel2 3 '(99 2)] '(1 1) )]                                                                                        ;
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

  [expect ( [neighborhood  1 2 '(100)] '#(0 1 2 3)  )
          ( [neighborhood 99 2 '(100)] '#(97 98 99) )
          ( [neighborhood  0 2 '(5 5)] '#(0 1 2 5 6 7 10 11 12) )
          ( [neighborhood 12 2 '(5 5)] [build-vector 25 id] )
          ( [neighborhood 24 2 '(5 5)] '#(12 13 14 17 18 19 22 23 24) )
          ( [neighborhood 50 2 '(1 100)] '#(48 49 50 51 52) )]
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
;; -- Compute --
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
  (define COL (lambda ps (build-vector (length ps) (lambda (sx) (make-synapse sx (list-ref ps sx))))))
  (define SP42 (default-sp '(4) '(2) '[global-inhibition . #t]
                `[potential-pools . ,(vector (COL 2001) (COL 0 2001 2001))]))
  [expect ( [connected-inputs SP42 0] #b0001 )  
          ( [connected-inputs SP42 1] #b0110 )]
                                                                                            ;
(define (boosted-overlaps sp overlaps)   ;; SP (ColVecOf OverlapX) -> (ColVecOf OverlapX)
  ;; produce overlaps adjusted by boost factors
  (vector-map
    (lambda (ov bf)
      (make-overlap (int<- (fx10k* bf (count ov))) (colx ov)))
    overlaps (sp-boost-factors sp)))
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

  [expect ((let ((SP5 (default-sp '(1) '(5) '[min-pct-overlap-duty-cycles . 100]
                        `[overlap-duty-cycles . ,'#(600 10000 30000 60000 5000)])))
             [min-duty-cycles-global SP5]) 600 )]  
                                                                                            ;
(define (min-duty-cycles-local sp cx)    ;; SP ColumnX -> DutyCycle
  ;; produce minodc of column as percent of max odc of columns in its neighborhood
  (let* ( (locality         (neighborhood cx (sp-inhibition-radius sp) (sp-column-dimensions sp)))
          (max-overlap-duty (fxvector-max (vector-refs (sp-overlap-duty-cycles sp) locality))))
    (fx10k* (sp-min-pct-overlap-duty-cycles sp) max-overlap-duty)))

  [expect ((let ((SP58 (default-sp '(5) '(8) '[inhibition-radius . 1] '[min-pct-overlap-duty-cycles . 2000]
                        `[overlap-duty-cycles . ,'#(7000 1000 5000  100 7800 5500 1000  10)])))
             [build-vector 8 (lambda (cx) [min-duty-cycles-local SP58 cx])]) 
                                               '#(1400 1400 1000 1560 1560 1560 1100 200) )]
                                                                                            ;
(define (update-duty-cycles-helper       ;; (ColVecOf DutyCycle) (ColVecOf Number) Nat -> (ColVecOf DutyCycle)
          duty-cycles new-input period)  
  ;; produce updated duty cycle estimate
  (vector-map (lambda (dc ni)
                (int<- (/ (+ (* dc (- period 1)) ni) period)))
              duty-cycles new-input))

  (let ((vec5 (lambda (x) (make-vector 5 x))))
    [expect ( [update-duty-cycles-helper (vec5 1000) (vec5    0) 1000] (vec5  999) )
            ( [update-duty-cycles-helper (vec5 1000) (vec5 1000) 1000] (vec5 1000) )
            ( [update-duty-cycles-helper (vec5 1000) '#(2000 4000 5000 6000 7000) 1000]
                                                     '#(1001 1003 1004 1005 1006)  )])
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

  [sp-overlap-duty-cycles-set! SP42 '#(0 0)]
  [sp-active-duty-cycles-set!  SP42 '#(0 0)]
  [expect ((begin [update-duty-cycles SP42 '#(0 0) '()] [sp-overlap-duty-cycles SP42]) '#(0 0) )]
                                                                                            ;
;; -- Inhibition radius --
                                                                                            ;
(define (avg-columns-per-input sp)       ;; SP -> Number
  ;; produce mean over dimensions of num-cols / num-inputs
  (list-average (map / (sp-column-dimensions sp) (sp-input-dimensions sp))))

  [expect ( [avg-columns-per-input SP42] 1/2 )]
                                                                                            ;
(define (avg-connected-span-1d sp cx)    ;; SP ColumnX -> Nat
  ;; produce width of range of input indices of connected synapses for column
  (let ((connected (vector-ref (sp-connected-synapses sp) cx)))
    (if (zero? connected)
        0
        (bitwise-span connected))))
  
  [expect ( [avg-connected-span-1d SP42 1] 2 )]
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
  
  (let ((SP2221 (default-sp '(2 2) '(2 1)
                  `[potential-pools . ,(vector (COL 0 2001 2001 0) (COL 0))])))
    [expect ( [avg-connected-span-2d SP2221 0] 3/2 )])
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

  [expect ( [inhibition-radius (default-sp '() '(57 31 2) '[global-inhibition . #t])] 57 )]
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

  [expect ((let ((SP84 (default-sp '(8) '(4) '[syn-perm-inactive-dec . 100] '[syn-perm-active-inc . 1000]
                   `[potential-pools . ,(vector (vector 02000 11200 20900 30400                        )
                                                (vector 01500                   41800 51200       74500)
                                                (vector             20140                   61100      )
                                                (vector 00400                               61780      ))])))
             [adapt-synapses SP84 #b01011001 '(0 1 2)]
             [sp-overlap-duty-cycles-set! SP84 '#(0 0 0 2500)] [sp-min-overlap-duty-cycles-set! SP84 '#(0 0 0 5000)]
             [bump-up-weak-columns SP84]
             [sp-potential-pools SP84]) (vector (vector 03000 11100 20800 31400                        )
                                                (vector 02500                   42800 51100       74400)
                                                (vector             20000                   62100      )
                                                (vector 00500                               61880      )))]
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

  [let ((SP11024 (default-sp '(1 1024) '(1 2048) '[potential-radius . 512])))
    (expect ( [map-column SP11024    0]    0 )
            ( [map-column SP11024 2047] 1023 ))]
  [let ((SP124 (default-sp '(12) '(4) '[potential-radius . 2]))
        (SP44 (default-sp '(4) '(4))))
    (do ((i 0 (add1 i))) ((= i 4)) (expect ( [map-column SP124 i] (add1 (* 3 i)) )))
    (do ((i 0 (add1 i))) ((= i 4)) (expect ( [map-column SP44  i] i )))]
  [let ((SP2D1 (default-sp '(36 12) '(12 4)))
        (SP2D2 (default-sp '(3 5) '(4 4))))
    (do ((i 0 (add1 i))) ((= i 5)) (expect ( [map-column SP2D1 (list-ref '( 0  4  5  7  47) i)]
                                                               (list-ref '(13 49 52 58 418) i))))
    (do ((i 0 (add1 i))) ((= i 3)) (expect ( [map-column SP2D2 (list-ref '(0 3 15) i)]
                                                               (list-ref '(0 4 14) i))))]
                                                                                            ;
(define (map-potential-v sp cx)          ;; SP ColumnX -> (vectorof InputX)
  ;; produce vector of potential input indices for column
  (let* ((center-input (map-column sp cx))
         (column-inputs (neighborhood center-input (sp-potential-radius sp) (sp-input-dimensions sp)))
         (num-potential (int<- (* (vector-length column-inputs) (sp-potential-pct sp)))))
    (vector-sample column-inputs num-potential)))

  (let ((sp (lambda (ni nc pr pp)
              (default-sp (list ni) (list nc) (cons 'potential-radius pr) (cons 'potential-pct pp)))))
    [expect ( (vector-sort < [map-potential-v (sp 12 4 2 1) 0])  '#(0 1 2 3)   )
            ( (vector-sort < [map-potential-v (sp 12 4 2 1) 2])  '#(5 6 7 8 9) )
            ( (vector-length [map-potential-v (sp 12 4 2 0.5) 0])  2           )
            ( [map-potential-v (sp 1 1 2 1.0) 0] '#(0) )])
  [let ((SP11024 (default-sp '(1 1024) '(1 2048) '[potential-radius . 2] '[potential-pct . 1])))
    (expect ( [map-potential-v SP11024 1023] '#(510 512 513 511 509) ))]
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
        
  (define SP16 (default-sp '(1) '(6) '[boost-strength . 10.0] 
                '[num-active-columns-per-inh-area . 1] '[inhibition-radius . 3]
                `[active-duty-cycles . ,(fx10k<- '#(0.1 0.3 0.02 0.04 0.7 0.12))]))
  [expect ([boost-factors-global SP16] '#(19477 2636 43348 35490 0048 15947) )]
                                                                                            ;
(define (boost-factors-local sp)         ;; SP -> (ColVecOf BoostFactor)
  ;; produce boost factors for all columns from local target density
  (vector-map-indexed
    (lambda (adc cx)                     ;; DutyCycle ColumnX -> BoostFactor
      (let* ( (mask-neighbors (neighborhood cx (sp-inhibition-radius sp) (sp-column-dimensions sp)))
              (target-density (/ (vector-average (vector-refs (sp-active-duty-cycles sp) mask-neighbors)) 10000.)))
        (boost-func sp adc target-density)))
    (sp-active-duty-cycles sp)))

  [expect ([boost-factors-local SP16] '#(11618 5066 69125 56595 97 27183) )]
                                                                                            ;
(define (calculate-overlap sp input-vec) ;; SP InputVec -> (ColVecOf OverlapX)
  ;; produce overlaps: vector of counts of intersection of input and synapses
  (build-vector 
    (sp-num-columns sp)
    (lambda (cx)
      (let ((count (bitwise-bit-count 
              (bitwise-and input-vec (vector-ref (sp-connected-synapses sp) cx)))))
        (make-overlap count cx)))))

  [expect ( [calculate-overlap SP42 #b1111] (vector (make-overlap 1 0) (make-overlap 2 1)) )
          ( [calculate-overlap SP42 #b1110] (vector (make-overlap 0 0) (make-overlap 2 1)) )
          ( [calculate-overlap SP42 #b0101] (vector (make-overlap 1 0) (make-overlap 1 1)) )]
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

  (let ((SP010 (default-sp '(0) '(10) '[stimulus-threshold . 6])))
    [expect ((list-sort < [inhibit-columns-global SP010 (list->overlaps '(1 2 1 4 8 3 12 5 4 1)) 0.3]) '(4 6)       )
            ((list-sort < [inhibit-columns-global SP010 (list->overlaps '(1 2 3 4 5 6 7 8 9 10)) 0.5]) '(5 6 7 8 9) )
            ((let ((l (list-sort < [inhibit-columns-global SP010 (list->overlaps '(1 1 1 7 7 7 7 8 9 10)) 0.5])))
                (and (member (car l) '(3 4 5)) (member (cadr l) '(4 5 6)) (equal? (cddr l) '(7 8 9)))) #t )])
                                                                                            ;
(define (tied-overlaps overlaps this-ov) ;; (vectorof OverlapX) (OverlapX) -> (vectorof ColumnX)
  ;; produce indices in neighborhood for overlaps with same count as this-ov
  (vector-filter id
    (vector-map-indexed
      (lambda (ov cx)
        (if (= (count ov) (count this-ov)) cx #f))
      overlaps)))

  (let ((mo (lambda (c x) (make-overlap c x))))
    [expect ( [tied-overlaps (list->overlaps '(1 2 3 1 3))       (mo 3 0)] '#(2 4) )
            ( [tied-overlaps (vector (mo 3 2) (mo 1 3) (mo 3 4)) (mo 3 0)] '#(0 2) )])
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
    
  (let ((SP010 (default-sp '(0) '(10) '[stimulus-threshold . 6])))
    [expect ((list-sort < [inhibit-columns-local SP010 (list->overlaps '(1 2 1 4 8 3 12 5 4 1)) 1.0]) '(4 6)       )
            ((list-sort < [inhibit-columns-local SP010 (list->overlaps '(1 2 3 4 5 6 7 8 9 10)) 1.0]) '(6 7 8 9) )
            ((list-sort < [inhibit-columns-local SP010 (list->overlaps '(1 1 1 7 7 7 7 8 9 10)) 0.5]) '() )])
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

;; ==================================== Examples ====================================

(define (random-bits size)
  (let ((w (expt 2 32)))
    (do ((i 1 (add1 i)) (r (random w) (+ (* w r) (random w)))) ((= i (div size 32)) r))))
                                                                                            ;
(define (percent-overlap x1 x2)
  (let ((min-x1x2 (min (bitwise-bit-count x1) (bitwise-bit-count x2))))
    (if (positive? min-x1x2)
        (* 100 (/ (bitwise-bit-count (bitwise-and x1 x2)) min-x1x2))
        0)))
                                                                                            ;
(define (percent->string x)
  (string-append (if (< x 10) "  " (if (< x 100) " " "")) 
                 (number->string (inexact (/ (int<- (* 10 x)) 10))) "%"))
                                                                                            ;
(define (corrupt-vector vector size noise-level)
  (let* ( (bit-xs (vector-sample (build-vector size id) (int<- (* noise-level size))))
          (flip   (vector->bitwise bit-xs)))
    (bitwise-xor vector flip)))
                                                                                            ;
(define (hello-sp) 
  (display "See nupic/examples/sp/hello_sp.py") (newline)
  (letrec* ( 
    (inp-dims '(32 32))
    (col-dims '(64 64))
    (ni   (apply * inp-dims))
    (sp   (make-sp* inp-dims col-dims '[global-inhibition . #t]
          `[potential-radius                . ,ni]
          `[num-active-columns-per-inh-area . ,(int<- (* 0.02 (apply * col-dims)))]
          `[syn-perm-active-inc             . ,(perm<- 0.01)]
          `[syn-perm-inactive-dec           . ,(perm<- 0.008)]))
    (create-input (lambda () (random-bits ni)))
    (run          (lambda (description input learn prev-cols)
                    (let* ( (cols  (list-sort < (compute sp input learn)))
                            (score (percent-overlap (list->bitwise cols) (list->bitwise prev-cols))))
                      (for-each display `(,description ": " ,(percent->string score) " ("))
                      (for-each (lambda (x) (display x) (display " ")) (take 15 cols)) 
                      (display "... ") (display (car (reverse cols))) (display ")") (newline)
                      cols)))
    (input1 (create-input))
    (input2 (create-input))
    (input3 (create-input))
    (input4 (corrupt-vector input3 ni 0.1))
    (input5 (corrupt-vector input4 ni 0.2)))
  (let* (
      (cols1 (run "Random 1" input1 #f '()))
      (cols2 (run "Random 2" input2 #f cols1))
      (cols3 (run "Random 3" input3 #t cols2))
      (cols4 (run "Repeat 3" input3 #t cols3)))
    (run "Noise .1" input4 #t cols4)
    (run "Noise .3" input5 #t cols4)))
  'ok)
                                                                                            ;
(define (sp-tutorial)
  (display "See nupic/examples/sp/sp_tutorial.py") (newline)
  (let* ( 
      (input-dimensions  '(1024 1))
      (column-dimensions '(2048 1))
      (input-size    (apply * input-dimensions))
      (column-number (apply * column-dimensions))
      (create-input  (lambda (_) (random-bits input-size)))
      (input-array   (create-input 0))
      (sp   (make-sp* input-dimensions column-dimensions
            `[potential-radius                . ,(int<- (* 0.5 input-size))]
            `[num-active-columns-per-inh-area . ,(int<- (* 0.02 column-number))]
            `[global-inhibition               . #t]
            `[syn-perm-active-inc             . ,(perm<- 0.01)]
            `[syn-perm-inactive-dec           . ,(perm<- 0.008)]))
      (active-cols   (compute sp input-array #f))
      (all-counts    (vector-map count (calculate-overlap sp input-array)))
      (active-counts (vector-refs all-counts (list->vector active-cols)))
      (mean          (lambda (vec) (int<- (vector-average vec))))
      (for-each-noise-level
        (lambda (description proc)
          (display description) (newline)
          (do ((i 0 (add1 i))) ((= i 11) )
            (let ((noise-level (/ i 10.0)))
              (for-each display
                `(,noise-level "  " ,(percent->string (proc noise-level)) #\newline)))))))
    (for-each display 
      `("Figure 1 - overlap count means" #\newline
        "all cols: " ,(mean all-counts) #\newline
        "  active: " ,(mean active-counts) #\newline))
    (for-each-noise-level "Figure 2 - noise:overlap linear" (lambda (nl)
      (percent-overlap input-array (corrupt-vector input-array input-size nl))))
    (for-each-noise-level "Figure 3 - without training" (lambda (nl)
      (percent-overlap (list->bitwise (compute sp input-array #f))
                       (list->bitwise (compute sp 
                         (corrupt-vector input-array input-size nl) #f)))))
    (let* ( (num-examples 10)
            (input-vectors (build-vector num-examples create-input))
            (epochs 30))
      (do ((epoch 0 (add1 epoch))) ((= epoch epochs))
        (do ((i 0 (add1 i))) ((= i num-examples))
          (compute sp (vector-ref input-vectors i) #t)))
      (for-each-noise-level "Figure 4 - with training: sigmoid" (lambda (nl)
        (percent-overlap (list->bitwise (compute sp (vector-ref input-vectors 0) #f))
                         (list->bitwise (compute sp 
                           (corrupt-vector (vector-ref input-vectors 0) input-size nl) #f))))))))

  #;(time (sp-tutorial))                 ;; uncomment to run on Chez load-program
