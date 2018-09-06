#!r6rs

;; === HTM-scheme/algorithms/htm_prelude Copyright 2017 Roger Turner. ===
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

(library (HTM-scheme HTM-scheme algorithms htm_prelude)
                                                                                            ;
(export
  add1
  make-list
  build-list
  take
  list-average
  unique
  int<-
  id
  build-vector 
  vector-filter
  vector-map-indexed
  vector-fold-left
  vector-average
  fxvector-max
  x10k
  fx10k<-
  fx10k*
  vector-extend
  vector-take
  vector-count
  vector-indices
  vector-refs
  list->bitwise
  vector->bitwise
  bitwise-span
  bitwise->list
  random
  random-seed!
  vector-sample
  key-word-args
  set-tracer!
  trace
  define-memoized)
                                                                                            ;
(import (rnrs))

;; -- Types --
;; Boolean, Number, Integer, Fixnum, (listof X), (vectorof X) ... = Scheme types
;; X, Y, Z      = type parameters (arbitrary types in function type specification etc)
;; X Y -> Z     = function with argument types X, Y and result type Z
;; Nat          = natural number (including zero) (Scheme Fixnum or exact Integer)
;; Fix10k       = Fixnum interpreted as number with 4 decimal places, ie 10000 = 1.0
;; Bits         = Integer (Bignum) interpreted bitwise

(define (add1 n)                         ;; Fixnum -> Fixnum
  (fx+ 1 n))
                                                                                            ;
(define (make-list n x)                  ;; Nat X -> (listof X)
  ;; produce list of n x's
  (vector->list (make-vector n x)))
                                                                                            ;
(define (build-list n proc)              ;; Nat (Nat -> X) -> (listof X)
  ;; produce list of n X's by applying proc to 0..n-1
  (let loop [ (i 0) (xs (list))]
    (cond [(fx=? i n) (reverse xs)]
          [else (loop (add1 i) (cons (proc i) xs))])))
                                                                                            ;
(define (take n xs)                      ;; Nat (listof X) -> (listof X)
  ;; produce first n elements of xs
  (let loop [ (n n) (xs xs) (ys (list))]
    (cond [(or (fxzero? n) (null? xs)) (reverse ys)]
          [else (loop (fx- n 1) (cdr xs) (cons (car xs) ys))])))
                                                                                            ;
(define (list-average l)                 ;; (listof Number) -> Number
  ;; produce mean of non-empty list
  (/ (apply + l) (length l)))
                                                                                            ;
(define (unique eql? xs)                 ;; (X X -> Boolean) (listof X) -> (listof X)
  ;; produce list with adjacent duplicates by eql? removed
  (cond [(null? xs) '()]
        [(null? (cdr xs)) xs]
        [(eql? (car xs) (cadr xs)) (unique eql? (cdr xs))]
        [else (cons (car xs) (unique eql? (cdr xs)))]))
                                                                                      ;
(define (int<- x)                        ;; Number -> Integer
  (exact (round x)))
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
                                                                                            ;
(define (vector-filter pred vec)         ;; (X -> Boolean) (vectorof X) -> (vectorof X)
  ;; produce vector of elements of vec for which (pred elt) is not false
  (list->vector (filter pred (vector->list vec))))
                                                                                            ;
(define (vector-map-indexed proc vec)    ;; (X Nat -> Y) (vectorof X) -> (vectorof Y)
  ;; produce vector by applying proc to each element of vec and its index
  (let ((result (make-vector (vector-length vec))))
    (do ((i 0 (add1 i))) ((fx=? i (vector-length vec)) result)
      (vector-set! result i (proc (vector-ref vec i) i)))))
                                                                                            ;
(define (vector-fold-left proc obj vec)  ;; (X Y -> X) X (vectorof Y) -> X
  ;; (proc ... (proc (proc obj vec-0) vec-1) ... vec-last)
  (let ((acc obj))
    (vector-for-each
      (lambda (x) (set! acc (proc acc x)))
      vec)
    acc))
                                                                                            ;
(define (vector-average vec)             ;; (vectorof Number) -> Number
  ;; (list-average (vector->list vec))
  (/ (vector-fold-left + 0 vec) (vector-length vec)))
                                                                                            ;
(define (fxvector-max vec)               ;; (vectorof Fixnum) -> Fixnum
  ;; produce max element of vec
  (vector-fold-left fxmax (least-fixnum) vec))
                                                                                            ;
(define x10k 10000)                      ;; Fixnum                                                                                            
                                                                                            ;
(define (fx10k<- x)                      ;; [vectorof] Number -> [vectorof] Fix10k
  ;; produce Fix10k representation of number or vector of numbers
  (if (vector? x) (vector-map fx10k<- x)
                  (int<- (* x x10k))))
                                                                                            ;
(define (fx10k* x y)                     ;; Fix10k Fix10k -> Fix10k
  ;; produce product
  (div (+ (* x y) (div x10k 2)) x10k))
                                                                                            ;
(define (vector-extend vec)              ;; (vectorof X) -> (vectorof X)
  ;; increase the size of vec
  (let* ( (len (vector-length vec))
          (result (make-vector (fxmax 1 (fx+ len (fxmin len 10000))))))
    (do ((i (fx- len 1) (fx- i 1))) ((fxnegative? i) result)
      (vector-set! result i (vector-ref vec i)))))
                                                                                            ;
(define (vector-take vec n)              ;; (vectorof X) Nat -> (vectorof X)
  ;; produce copy of first n elements of vec, or copy of vec if n > length
  (let* ( (size (fxmin n (vector-length vec)))
          (result (make-vector size)))
    (do ((i 0 (add1 i))) ((fx=? i size) result)
      (vector-set! result i (vector-ref vec i)))))
                                                                                            ;
(define (vector-count proc vec)          ;; (X -> Boolean) (vectorof X) -> Nat
  ;; produce count of elements of vec for which (proc elt) is not false
  (vector-fold-left 
    (lambda (acc elt) (if (proc elt) (add1 acc) acc))
    0 vec))
                                                                                            ;
(define (vector-indices vec)             ;; (vectorof X) -> (listof Nat)
  ;; produce list of indices for which vector element is not false
  (do [ (i 0 (add1 i)) 
        (l (list) (if (vector-ref vec i) (cons i l) l))]
      ((fx=? i (vector-length vec)) l)))
                                                                                            ;
(define (vector-refs vec refs)           ;; (vectorof X) (vectorof Nat) -> (vectorof X)
  ;; produce vector of selected elements of vec indexed by refs
  (let ((vrefs (make-vector (vector-length refs))))
    (do ((i 0 (add1 i))) ((fx=? i (vector-length refs)) vrefs)
      (vector-set! vrefs i (vector-ref vec (vector-ref refs i))))))
                                                                                            ;
(define (list->bitwise ns)               ;; (listof Nat) -> Bits
  ;; produce bitwise value by setting bits indexed by elements of ns
  (fold-left (lambda (acc n) (bitwise-copy-bit acc n 1)) 0 ns))
                                                                                            ;
(define (vector->bitwise vec)            ;; (vectorof Nat) -> Bits
  (list->bitwise (vector->list vec)))
                                                                                            ;
(define (bitwise-span bits)              ;; Bits -> Nat
  ;; produce span of 1 bits in bits
  (fx- (bitwise-length bits) (bitwise-first-bit-set bits)))
                                                                                            ;
(define (bitwise->list bits)             ;; Bits -> (listof Nat)
  ;; produce indices in ascending order of set bits in bits
  ;; this is much faster in Racket 7.0 than loop below
  (let loop ((index (bitwise-length bits)) (result (list)))
    (if (fxnegative? index) result
      (loop (fx- index 1)
            (if (bitwise-bit-set? bits index)
              (cons index result)
              result))))
  ;; this is faster in Chez 9.5 than loop above
  #;(let loop ((bits bits) (result (list)))
    (if (positive? bits)
      (let ((b (bitwise-first-bit-set bits)))
        (loop (bitwise-copy-bit bits b 0) (cons b result)))
      (reverse result)))
  )
                                                                                            ;
(define (random n)                       ;; Nat -> Nat
  ;; produce random integer in range 0..n-1
  ;; alternative?: (set! *random-state* (mod (+ 1013904223 (* *random-state* 1664525)) 4294967296))
  (let ((Lehmer-modulus   2147483647)
        (Lehmer-multiplier     48271))
    (set! *random-state* (mod (* *random-state* Lehmer-multiplier) Lehmer-modulus))
    (mod *random-state* n)))
    
  (define *random-state* 48271)
                                                                                            ;
(define (random-seed! n)                 ;; Nat ->
  (set! *random-state* n))
                                                                                            ;
(define (vector-sample source size)      ;; (vectorof X) Nat -> (vectorof X)
  ;; produce random selection of length size from source (using Durstenfeld shuffle)
  (let* ( (source-length (vector-length source))
          (shuffle (if (fx<? size source-length) size (fx- size 1))))
    (cond [(fxzero? size) '#()]
          [(fx=? size 1)
           (vector (vector-ref source (random source-length)))]
          [else
            (let ((vec (vector-take source source-length)))
              (do ((n 0 (add1 n))) ((fx=? n shuffle) (vector-take vec size))
                (let* ((r (fx+ n (random (fx- source-length n))))
                       (t (vector-ref vec r)))
                  (vector-set! vec r (vector-ref vec n))
                  (vector-set! vec n t))))])))
                                                                                            ;
(define (key-word-args args defaults check-args)    ;; (listof KWarg) (listof KWarg) -> (listof X)
  ;; KWarg is [key . value]; produce list of default values overridden by arg values
  (when check-args
    (for-each (lambda (arg-kv)
                (unless (assq (car arg-kv) defaults)
                  (error #f "unknown keyword arg" arg-kv)))
              args))
  (map (lambda (default-kv)
         (let ((kv (assq (car default-kv) args)))  ;; if this default in args
           (if kv (cdr kv) (cdr default-kv))))     ;; then use given val
       defaults))
  
(define *tracer* (lambda x (if #f #f)))

(define (set-tracer! proc)
  (set! *tracer* proc))
                                                                                         ;
(define trace
  (lambda x (apply *tracer* x)))
  
(define-syntax define-memoized           ;; Function-defn -> defn with arg/result cache
  (syntax-rules ()
    ((define-memoized (f arg ...) body ...)
      (define f
        (let ((cache (list)))
          (lambda (arg ...)
            (cond ((assoc `(,arg ...) cache) => cdr)
                  (else (let ((val (begin body ...)))
                          (set! cache (cons (cons `(,arg ...) val) cache))
                          val)))))))))

)