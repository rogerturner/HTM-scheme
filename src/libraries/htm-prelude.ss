#!r6rs
;; ============ HTM-scheme Prelude Copyright 2017 Roger Turner. ============
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

(library (libraries htm-prelude)

  (export
    expect
    list-average
    add1
    build-list
    take
    make-list
    list->bitwise
    int<-
    id
    build-vector 
    id
    vector-filter
    vector-map-indexed
    vector-fold-left
    vector-average
    fxvector-max
    x10k
    fx10k<-
    fx10k*
    vector-take
    vector-count
    vector-indices
    vector-refs
    vector->bitwise
    bitwise-span
    random 
    vector-sample
    key-word-args)
    
  (import (rnrs))

(define-syntax expect                    ;; ((X ... -> Y) X ...) Y -> [error]
  ;; check that function application(s) to arguments match expected values
  (lambda (x)                            
    (syntax-case x ()                                  ;; [expect ([fn args] expected ) ... ]
      [ (_ (expr expected) ...)                        ;; expr matches [fn args]
        #'(define (expect) #f)                         ;; disable tests in library
  #;    #'(begin (let ((result expr))                  ;; eval expr just once
                   (unless (equal? result expected)    ;; no output if check passes
                           (for-each display `("**" expr #\newline 
                             "  expected: " expected #\newline 
                             "  returned: " ,result  #\newline)))) ...)])))
                                                                                            ;
(define (list-average l)                 ;; (listof Number) -> Number
  ;; produce mean of non-empty list
  (/ (apply + l) (length l)))

  [expect ( [list-average (list 1 2 27)]  10 )]
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
                (let* ((r (random (fx- source-length n)))
                       (t (vector-ref vec r)))
                  (vector-set! vec r (vector-ref vec n))
                  (vector-set! vec n t))))])))

  [define SOURCE [build-vector 100 id]]
  [define SELECT [vector->list [vector-sample SOURCE 50] ]]
  [expect ( (length SELECT) 50 )
          ( (list? (for-all memq SELECT [make-list 50 [vector->list SOURCE]])) #t )  ;; all from source
          ( (map (lambda (x) (length (remq x SELECT))) SELECT) [make-list 50 49]  )] ;; all different
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

)