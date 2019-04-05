#!chezscheme

;; === HTM-scheme/algorithms/htm_prelude Copyright 2019 Roger Turner. ===
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

  ;; Library functions specialized for HTM-scheme
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.

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
  indexes
  vector-map-x
  vector-for-each-x
  vector-fold-left
  vector-average
  fxvector-max
  fx3
  fx3<-
  fx3*
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
  fxsearch
  key-word-args
  intersect1d
  union1d
  setdiff1d
  in1d
  include-by-mask
  exclude-by-mask
  define-memoized)
                                                                                            ;
(import
  (except (chezscheme) add1 make-list random reset))

;; -- Types --
;; Boolean, Number, Integer, Fixnum, (Listof X), (Vectorof X) ... = Scheme types
;; X, Y, Z  = type parameters (arbitrary types in function type specification etc)
;; X Y -> Z = function with argument types X, Y and result type Z; [X] = optional X
;; Nat      = natural number (including zero) (Scheme Fixnum or exact Integer)
;; Fixnum3  = Fixnum interpreted as number with 3 decimal places, ie 1000 = 1.0
;; Bits     = Integer (Bignum) interpreted bitwise

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
    (cond [(fx=? i n) (reverse! xs)]
          [else (loop (add1 i) (cons (proc i) xs))])))
                                                                                            ;
(define (take n xs)                      ;; Nat (listof X) -> (listof X)
  ;; produce first n elements of xs
  (let loop [ (n n) (xs xs) (ys (list))]
    (cond [(or (fxzero? n) (null? xs)) (reverse! ys)]
          [else (loop (fx- n 1) (cdr xs) (cons (car xs) ys))])))
                                                                                            ;
(define (list-average l)                 ;; (listof Number) -> Number
  ;; produce mean of non-empty list
  (/ (apply + l) (length l)))
                                                                                            ;
(define (unique eql? xs)                 ;; (X X -> Boolean) (listof X) -> (listof X)
  ;; produce list with adjacent duplicates by eql? removed
  ;;  (cond [(null? xs) '()]
  ;;        [(null? (cdr xs)) xs]
  ;;        [(eql? (car xs) (cadr xs)) (unique eql? (cdr xs))]
  ;;        [else (cons (car xs) (unique eql? (cdr xs)))]))
  (if (null? xs)  xs
    (let loop ((first xs) (next (cdr xs)))
      (cond
        [(null? next) xs]
        [(eql? (car first) (car next))
          (let skip ((next (cdr next)))
            (cond
              [(null? next)
                (set-cdr! first next)
                xs]
              [(eql? (car first) (car next))
                (skip (cdr next))]
              [else 
                (set-cdr! first next)
                (loop next (cdr next))]))]
        [else
          (loop next (cdr next))]))))
                                                                                      ;
(define (int<- x)                        ;; Number -> Integer
  (exact (round x)))
                                                                                            ;
(define (id x)                           ;; X -> X
  x)
                                                                                            ;
(define (build-vector n f)               ;; Nat (Nat -> X) -> (vectorof X)
  ;; produce vector length n by applying f to indices
  (let ((v (make-vector n)))
    (do ((i 0 (add1 i))) ((fx=? i n) v)
      (vector-set! v i (f i)))))
                                                                                            ;
(define (vector-filter pred vec)         ;; (X -> Boolean) (vectorof X) -> (vectorof X)
  ;; produce vector of elements of vec for which (pred elt) is not false
  (list->vector (filter pred (vector->list vec))))
                                                                                            ;
(define (indexes seq)                    ;; (vectorof X) -> (vectorof Nat)                                                                 
  ;; produce indexes of vector or list   ;; (listof X)   -> (listof Nat)
  (if (vector? seq)
      (build-vector (vector-length seq) id)
      (build-list   (length seq)        id)))
                                                                                            ;
(define (with-index folder f vs)         ;; ((? -> ?) ?... -> ?) (?... -> ?) (listof (vectorof ?)) -> ?
  (apply folder
    (lambda xs (apply f xs))
    (append vs (list (indexes (car vs))))))
                                                                                            ;
(define (vector-map-x f . vs)            ;; (X ... Nat -> Y) (vectorof X) ... -> (vectorof Y)
  ;; produce vector by applying f to each element of vs and its index
  (with-index vector-map f vs))
                                                                                            ;
(define (vector-for-each-x f . vs)       ;; (X ... Nat -> ) (vectorof X) ... -> 
  ;; apply f to each element of vs and its index
  (with-index vector-for-each f vs))
                                                                                            ;
(define (vector-fold-left f o . vs)      ;; (X Y ... -> X) X (vectorof Y) ... -> X
  (let ((acc o))
    (apply vector-for-each
      (lambda xs (set! acc (apply f acc xs)))
      vs)
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
(define fx3 1000)                        ;; Fixnum                                                                                            
                                                                                            ;
(define (fx3<- x)                        ;; [Vectorof] Number -> [Vectorof] Fixnum3
  ;; produce Fixnum3 representation of number or vector of numbers
  (if (vector? x) (vector-map fx3<- x)
                  (int<- (* x fx3))))
                                                                                            ;
(define (fx3* x y)                       ;; Fixnum3 Fixnum3 -> Fixnum3
  ;; produce product
  (div (+ (* x y) (div fx3 2)) fx3))
                                                                                            ;
(define (vector-extend vec)              ;; (vectorof X) -> (vectorof X)
  ;; increase the size of vec
  (let* ( (len (vector-length vec))
          (result (make-vector (fxmax 1 (fx+ len (fxmin len 10000))))))
    (do ((i (fx- len 1) (fx- i 1))) ((fxnegative? i) result)
      (vector-set! result i (vector-ref vec i)))))
                                                                                            ;
(define (vector-take n vec)              ;; Nat (vectorof X) -> (vectorof X)
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
  #;(let loop ((index (bitwise-length bits)) (result (list)))
    (if (fxnegative? index) result
      (loop (fx- index 1)
            (if (bitwise-bit-set? bits index)
              (cons index result)
              result))))
  ;; this is faster in Chez 9.5 than loop above
  ;; negative bits (infinite 1s) is clipped to fixnum range!
  (let loop (
        (bits (if (negative? bits)
                  (bitwise-and bits (greatest-fixnum))
                  bits))
        (result (list)))
    (if (zero? bits)  (reverse! result)
      (let ((b (bitwise-first-bit-set bits)))
        (loop (bitwise-copy-bit bits b 0) (cons b result)))))
  )
                                                                                            ;
(define (random n)                       ;; Nat -> Nat
  ;; produce random integer in range 0..n-1
  ;; alternative?: (set! *random-state* (mod (+ 1013904223 (* *random-state* 1664525)) 4294967296))
  (let ((Lehmer-modulus   2147483647)
        (Lehmer-multiplier     48271))
    (set! *random-state* (mod (* *random-state* Lehmer-multiplier) Lehmer-modulus))
    (mod *random-state* n)))
                                                                                            ;
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
            (let ((vec (vector-take source-length source)))
              (do ((n 0 (add1 n))) ((fx=? n shuffle) (vector-take size vec))
                (let* ((r (fx+ n (random (fx- source-length n))))
                       (t (vector-ref vec r)))
                  (vector-set! vec r (vector-ref vec n))
                  (vector-set! vec n t))))])))
                                                                                            ;
(define (fxsearch v target)              ;; (vectorof Fixnum) Fixnum -> Fixnum|#f
  (let search ((left 0) (right (fx- (vector-length v) 1)))
    (if (fx>? left right) #f
      (let* ( (mid (fxdiv (fx+ left right) 2))
              (element (vector-ref v mid))) 
        (cond 
          [ (fx<? element target) (search (add1 mid) right) ]
          [ (fx<? target element) (search left (fx- mid 1)) ]
          [else element])))))
                                                                                            ;
(define (key-word-args args defaults)    ;; (listof KWarg) (listof KWarg) -> (listof X)
  ;; KWarg is [key . value]; produce list of default values overridden by arg values
  (map (lambda (default-kv)
         (let ((kv (assq (car default-kv) args)))  ;; if this default in args
           (if kv (cdr kv) (cdr default-kv))))     ;; then use given val
       defaults))
                                                                                            ;
(define (intersect1d l1 l2)              ;; {Fixnum} {Fixnum} -> {Fixnum}
  ;; produce intersection of sorted lists of fixnums
  ;; (assert (equal? l1 (list-sort fx<? l1)))
  ;; (assert (equal? l2 (list-sort fx<? l2)))
  (let loop ((l1 (unique fx=? l1)) (l2 (unique fx=? l2)) (result (list)))
    (cond [(or (null? l1) (null? l2)) (reverse! result)]
          [(fx<? (car l1) (car l2)) (loop (cdr l1) l2 result)]
          [(fx<? (car l2) (car l1)) (loop l1 (cdr l2) result)]
          [else (loop (cdr l1) (cdr l2) (cons (car l1) result))])))
                                                                                            ;
(define (union1d l1 l2)                  ;; {Fixnum} {Fixnum} -> {Fixnum}
  ;; produce union of sorted lists of fixnums
  ;; (assert (equal? l1 (list-sort fx<? l1)))
  ;; (assert (equal? l2 (list-sort fx<? l2)))
  (let loop ((l1 l1) (l2 l2) (result (list)))
    (cond [(and (null? l1) (null? l2)) (unique fx=? (reverse! result))]
          [(null? l1) (append (reverse! result) l2)]
          [(null? l2) (append (reverse! result) l1)]
          [(fx<? (car l1) (car l2)) (loop (cdr l1) l2 (cons (car l1) result))]
          [(fx<? (car l2) (car l1)) (loop l1 (cdr l2) (cons (car l2) result))]
          [else (loop (cdr l1) (cdr l2) (cons (car l1) result))])))
                                                                                            ;
(define (setdiff1d l1 l2)                ;; {Fixnum} {Fixnum} -> {Fixnum}
  ;; produce difference of sorted lists of fixnums
  ;; (assert (equal? l1 (list-sort fx<? l1)))
  ;; (assert (equal? l2 (list-sort fx<? l2)))
  (let loop ((l1 (unique fx=? l1)) (l2 (unique fx=? l2)) (result (list)))
    (cond [(null? l1) (reverse! result)]
          [(null? l2) (append (reverse! result) l1)]
          [(fx<? (car l1) (car l2)) (loop (cdr l1) l2 (cons (car l1) result))]
          [(fx<? (car l2) (car l1)) (loop l1 (cdr l2) (cons (car l1) result))]
          [else (loop (cdr l1) (cdr l2) result)])))
                                                                                            ;
(define (in1d x1s x2s)                   ;; {X} {X} -> Bits
  ;; produce index mask of x1s that are also in x2s
  (fold-left 
    (lambda (acc e1 e1x)
      (if (memv e1 x2s)
        (bitwise-copy-bit acc e1x 1)
        acc))
    0 x1s (build-list (length x1s) id)))
                                                                                            ;
(define (include-by-mask xs mask)        ;; {X} Bits -> {X}
  ;; extract elements of xs corresponding to 1 bits in mask
  (vector->list
    (vector-refs 
      (list->vector xs) 
      (list->vector (bitwise->list mask)))))
                                                                                            ;
(define (exclude-by-mask xs mask)        ;; {X} Bits -> {X}
  (vector->list
    (vector-refs 
      (list->vector xs) 
      (list->vector 
        (bitwise->list
          (bitwise-bit-field (bitwise-not mask) 0 (length xs)))))))
                                                                                            ;
  ;; List comprehensions from the Programming Praxis standard prelude
  ;; Clauses: (see https://programmingpraxis.com/contents/standard-prelude/ for details)
  ;;   (var range [first] past [step]) — Bind var to first, first + step, ...
  ;;   (var in list)                   — Loop over elements of list
  ;;   (var is expr)                   — Bind var to value of expr
  ;;   (pred? expr)                    — Include elements for which (pred? x) is non-#f.
                                                                                            ;
(define-syntax fold-of                   ;; Op Base Expr Clause ...
  (syntax-rules (range in is)
    ((_ "z" f b e) (set! b (f b e)))
    ((_ "z" f b e (v range fst pst stp) c ...)
      (let* ((x fst) (p pst) (s stp)
             (le? (if (positive? s) <= >=)))
        (do ((v x (+ v s))) ((le? p v) b)
          (fold-of "z" f b e c ...))))
    ((_ "z" f b e (v range fst pst) c ...)
      (let* ((x fst) (p pst) (s (if (< x p) 1 -1)))
        (fold-of "z" f b e (v range x p s) c ...)))
    ((_ "z" f b e (v range pst) c ...)
      (fold-of "z" f b e (v range 0 pst) c ...))
    ((_ "z" f b e (x in xs) c ...)
      (do ((t xs (cdr t))) ((null? t) b)
        (let ((x (car t)))
          (fold-of "z" f b e c ...))))
    ((_ "z" f b e (x is y) c ...)
      (let ((x y)) (fold-of "z" f b e c ...)))
    ((_ "z" f b e p? c ...)
      (if p? (fold-of "z" f b e c ...)))
    ((_ f i e c ...)
      (let ((b i)) (fold-of "z" f b e c ...)))))
                                                                                            ;
(define-syntax list-of (syntax-rules ()  ;; Expr Clause ...
  ((_ arg ...) (reverse (fold-of
    (lambda (d a) (cons a d)) '() arg ...)))))
                                                                                            ;
(define-syntax sum-of (syntax-rules ()   ;; Expr Clause ...
  ((_ arg ...) (fold-of + 0 arg ...))))
                                                                                            ;
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

;; -- Smoke tests --
                                                                                            ;
(define-syntax expect                    ;; ((X ... -> Y) X ...) Y -> [error]
  ;; check that function application(s) to arguments match expected values
  (lambda (x)                            
    (syntax-case x ()                    ;; [expect ([fn args] expected ) ... ]
      [ (_ (expr expected) ...)          ;; expr matches [fn args]
        #'(begin 
            (let ((result expr))         ;; eval expr just once, no output if check passes
              (unless (equal? result expected)
                (for-each display 
                  `("**" expr #\newline 
                    "  expected: " expected #\newline 
                    "  returned: " ,result  #\newline))
                (exit))) ...)])))
                                                                                            ;
  [expect ( [list-average (list 1 2 27)]  10 )]
  ;[expect ( [unique fx=? '(1 2 2 3 4 4)] '(1 2 3 4) )]
  [expect ( [build-list 3 -] '(0 -1 -2) )]
  [expect ( [take 2 '(a 2 c)] '(a 2) )]
  [expect ( [make-list 0 9] '()      )
          ( [make-list 3 9] '(9 9 9) )]
  [expect ( [list->bitwise '(0 2 3)] #b1101 )]
  [expect ( [int<- 2.500000000000001] 3 )
          ( [int<- 5/2              ] 2 )
          ( [int<- 7/2              ] 4 )]
  [expect ( [build-vector 3 id] '#(0 1 2) )]
  [expect ( [vector-filter even? '#(1 2 3 4 5)] '#(2 4) )]
  [expect ( [indexes '(1 2 3)]  '(0 1 2)  )]
  [expect ( [indexes '#(1 2 3)] '#(0 1 2) )]
  [expect ( [vector-map-x (lambda (x y i) (+ x y i)) '#(10 11 12) '#(1 1 1)] '#(11 13 15) )]
  [expect ( [vector-fold-left (lambda (l x) (cons x l)) '() '#(0 1 2)] '(2 1 0) )]
  [expect ( [vector-fold-left (lambda (s x y) (* s (+ x y))) 1 '#(0 1 2) '#(1 2 3)] 15 )]
  [expect ( [vector-average '#(1 2 5)] 8/3 )]
  [expect ( [fxvector-max '#(0 -3 7 5)] 7 )]
  [expect ( [fx3<- '#(1/3 1.0)] '#(333 1000) )]
  [expect ( [fx3* 500 500] 250 )]
  [expect ( [vector-take 0 '#(0 1 2 3)] '#()       )
          ( [vector-take 2 '#(0 1 2 3)] '#(0 1)    )
          ( [vector-take 5 '#(0 1 2 3)] '#(0 1 2 3))]
  [expect ( [vector-count zero? '#(0 1 2 0 3)] 2 )
          ( [vector-count not   '#(#t 0 add1)] 0 )]
  [expect ( [vector-indices '#(#f #f #f)] '()    )
          ( [vector-indices '#(#t #f  0)] '(2 0) )]
  [expect ( [vector-refs '#(0 11 22 33 44) '#(0 2 4)] '#(0 22 44) )]
  [expect ( [bitwise-span #b01011101000] 7 )]
  [let* ( (source (build-vector 100 id))
          (select (vector->list [vector-sample source 50] )))
    (expect ( (length select) 50 )
            ( (list? (for-all memq select (make-list 50 (vector->list source)))) #t )  ;; all from source
            ( (map (lambda (x) (length (remq x select))) select) (make-list 50 49)) )] ;; all different
  [expect ( [(lambda (x . args)          ;; example fn with kws k1 & k2: default is (x 1 2)
              (append (list x) [key-word-args args '([k1 . 1] [k2 . 2])] ))  ;; end of fn defn
            99 '[k2 . 22] '[k1 . 11] #;'[k3 . 33] ] '(99 11 22) )]  ;; apply fn to 99 '[k2 ...
  [expect ( [intersect1d     '(1 2 3 4) '(1 3 5)] '(1 3) )]
  [expect ( [setdiff1d       '(1 2 3 4) '(1 3 5)] '(2 4) )]
  [expect ( [in1d            '(1 2 3 4) '(1 3 5)]  #b101 )]
  [expect ( [include-by-mask '(1 2 3 4)  #b1101]  '(1 3 4))]
  [expect ( [exclude-by-mask '(1 2 3 4)  #b1101]  '(2))]
  
)