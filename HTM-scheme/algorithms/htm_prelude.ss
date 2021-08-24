;; === HTM-scheme Prelude  (C) 2019-2021 Roger Turner. ===
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
  #|
  
Utility functions specialized for HTM-scheme (Fixnums used wherever possible)
Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.

"Lisp's parentheses are the bumps on the top of Lego" [Paul Graham]
"The only honest function names are car and cdr"  [citation needed]

  |#
  
  #!chezscheme

(library (HTM-scheme HTM-scheme algorithms htm_prelude)
                                                                                            ;
(export
  build-list
  take
  list-average
  list-refs!
  unique!
  condense!
  argmax-multi
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
  vector-sample
  u32-sample
  fxsearch
  key-word-args
  do-with-progress
  intersect1d
  union1d
  setdiff1d
  in1d
  include-by-mask
  exclude-by-mask
  thread-limit
  threaded-vector-for-each
  cost-center-1
  cost-center-2
  define-memoized)
                                                                                            ;
(import (chezscheme))
                                                                                            ;
  (implicit-exports #f)

;; -- Types --
;; Boolean, Number, Integer, Fixnum, (Listof X), (Vectorof X) ... = Scheme types
;; X, Y, Z  = type parameters (arbitrary types in function type specification etc)
;; X Y -> Z = function with argument types X, Y and result type Z; [Y] = optional Y
;; {X}      = abbreviation for (Listof X)
;; Nat      = natural number (including zero) (Scheme Fixnum or exact Integer)
;; Fixnum3  = Fixnum interpreted as number with 3 decimal places, ie 1000 = 1.0
;; Bits     = Integer (Bignum) interpreted bitwise
;; KWarg    = Pair (key . value)

(define (build-list n proc)              ;; Nat (Nat -> X) -> (Listof X)
  ;; produce list of n X's by applying proc to 0..n-1 (in reverse order)
  (let build ([i (fx1- n)] [l (list)])
    (if (fxnegative? i) l
      (build (fx1- i) (cons (proc i) l)))))
                                                                                            ;
(define (take n xs)                      ;; Nat (Listof X) -> (Listof X)
  ;; produce first n elements of xs
  (let recur ([n n] [xs xs])
    (if (fxzero? n)  '()
      	(cons (car xs)
      	      (recur (fx1- n) (cdr xs))))))
                                                                                            ;
(define (list-average l)                 ;; (Listof Number) -> Number
  ;; produce mean of non-empty list
  (/ (apply + l) (length l)))
                                                                                            ;
(define (list-refs! ls xs)               ;; {X} {Nat} -> {X}
  ;; mutate ls by removing elements not in xs, which are sorted indexes of ls elements
  (if (null? xs)  '()
    (let ([ls (list-tail ls (car xs))])
      (let next-x ([prev ls] [xs (cdr xs)] [prev-x (car xs)])
        (if (null? xs)  (begin (set-cdr! prev '())  ls)
          (let ([this-x (car xs)])
            (unless (fx=? this-x (fx1+ prev-x))
              (set-cdr! prev (list-tail prev (fx- this-x prev-x))))
            (next-x (cdr prev) (cdr xs) this-x)))))))
                                                                                            ;
(define (unique! eql? xs)                ;; (X X -> Boolean) (Listof X) -> (Listof X)
  ;; mutate xs removing adjacent duplicates by eql?
  ;;  (cond [(null? xs) '()]
  ;;        [(null? (cdr xs)) xs]
  ;;        [(eql? (car xs) (cadr xs)) (unique eql? (cdr xs))]
  ;;        [else (cons (car xs) (unique eql? (cdr xs)))]))
  (if (null? xs)  xs
    (let loop ([first xs] [next (cdr xs)])
      (cond
        [(null? next) xs]
        [(eql? (car first) (car next))
          (let skip ([next (cdr next)])
            (cond
              [(null? next)
                (set-cdr! first '())
                xs]
              [(eql? (car first) (car next))
                (skip (cdr next))]
              [else 
                (set-cdr! first next)
                (loop next (cdr next))]))]
        [else
          (loop next (cdr next))]))))
                                                                                            ;
(define (condense! l)                    ;; (Listof X) -> (Listof X)
  ;; mutate l by removing eq? duplicates
  (let condense ([len (length l)] [l l])
    (if (fx<? len 2)  l
      (let* ( [half-len (fxdiv len 2)]
              [last (list-tail l (fx1- half-len))]
              [rest (cdr last)])
        (set-cdr! last '())
        (let* ( [filtered-begin (condense half-len l)]
                [filtered-end
                  (let ([l (filter (lambda (x)
                               (not (memq x filtered-begin)))
                             rest)])
                    (condense (length l) l))])
          (append! filtered-begin filtered-end))))))
                                                                                            ;
(define (argmax-multi a group-keys)      ;; {Fixnum} {Fixnum} -> {Nat}
  ;; Get the indices of the max values of each group in 'a', grouping the
  ;; elements by their corresponding value in group-keys
  (let per-group ([a a] [group-keys group-keys] [index 0] [result (list)])
    (cond
      [ (null? group-keys) (reverse! result) ]
      [ else
        (let ([group (car group-keys)])
          (let next ([a a] [group-keys group-keys] [index index]
                     [max-in-group (least-fixnum)] [index-of-max #f])
            (cond
              [ (null? group-keys) (per-group a group-keys index (cons index-of-max result))]
              [ (fx=? (car group-keys) group)
                  (if (fx>? (car a) max-in-group)
                    (next (cdr a) (cdr group-keys) (fx1+ index) (car a) index)
                    (next (cdr a) (cdr group-keys) (fx1+ index) max-in-group index-of-max)) ]
              [ else (per-group a group-keys index (cons index-of-max result)) ]))) ])))
                                                                                            ;
(define (int<- x)                        ;; Number -> Integer
  (exact (round x)))
                                                                                            ;
(define (id x)                           ;; X -> X
  x)
                                                                                            ;
(define (build-vector n f)               ;; Nat (Nat -> X) -> (Vectorof X)
  ;; produce vector length n by applying f to indices
  (let ([v (make-vector n)])
    (do ([i 0 (fx1+ i)]) ((fx=? i n) v)
      (vector-set! v i (f i)))))
                                                                                            ;
(define (vector-filter pred vec)         ;; (X -> Boolean) (Vectorof X) -> (Vectorof X)
  ;; produce vector of elements of vec for which (pred elt) is not false
  (list->vector (filter pred (vector->list vec))))
                                                                                            ;
(define (indexes seq)                    ;; (Vectorof X) -> (Vectorof Nat)                                                                 
  ;; produce indexes of vector or list   ;; (Listof X)   -> (Listof Nat)
  (if (vector? seq)
    (build-vector (vector-length seq) id)
    (build-list   (length seq)        id)))
                                                                                            ;
(define (with-index folder f vs)         ;; ((? -> ?) ?... -> ?) (?... -> ?) (Listof (Vectorof ?)) -> ?
  (apply folder
    (lambda xs (apply f xs))
    (append vs (list (indexes (car vs))))))
                                                                                            ;
(define (vector-map-x f . vs)            ;; (X ... Nat -> Y) (Vectorof X) ... -> (Vectorof Y)
  ;; produce vector by applying f to each element of vs and its index
  (with-index vector-map f vs))
                                                                                            ;
(define (vector-for-each-x f . vs)       ;; (X ... Nat -> ) (Vectorof X) ... -> 
  ;; apply f to each element of vs and its index
  (with-index vector-for-each f vs))
                                                                                            ;
(define (vector-fold-left f o . vs)      ;; (X Y ... -> X) X (Vectorof Y) ... -> X
  (let ([acc o])
    (apply vector-for-each
      (lambda xs (set! acc (apply f acc xs)))
      vs)
    acc))
                                                                                            ;
(define (vector-average vec)             ;; (Vectorof Number) -> Number
  ;; (list-average (vector->list vec))
  (/ (vector-fold-left + 0 vec) (vector-length vec)))
                                                                                            ;
(define (fxvector-max vec)               ;; (Vectorof Fixnum) -> Fixnum
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
  (fxdiv (fx+ (fx* x y) (fxdiv fx3 2)) fx3))
                                                                                            ;
(define (vector-extend vec)              ;; (Vectorof X) -> (Vectorof X)
  ;; increase the size of vec
  (let* ( [len    (vector-length vec)]
          [result (make-vector (fxmax 2 (fx+ len (fxdiv len 2))))])
    (do ([i (fx1- len) (fx1- i)]) ((fxnegative? i) result)
      (vector-set! result i (vector-ref vec i)))))
                                                                                            ;
(define (vector-take n vec)              ;; Nat (Vectorof X) -> (Vectorof X)
  ;; produce copy of first n elements of vec, or copy of vec if n > length
  (let* ( [size   (fxmin n (vector-length vec))]
          [result (make-vector size)])
    (do ([i 0 (fx1+ i)]) ((fx=? i size) result)
      (vector-set! result i (vector-ref vec i)))))
                                                                                            ;
(define (vector-count proc vec)          ;; (X -> Boolean) (Vectorof X) -> Nat
  ;; produce count of elements of vec for which (proc elt) is not false
  (vector-fold-left 
    (lambda (acc elt) (if (proc elt) (fx1+ acc) acc))
    0 vec))
                                                                                            ;
(define (vector-indices vec)             ;; (Vectorof X) -> (Listof Nat)
  ;; produce list of indices for which vector element is not false
  (do ( [i 0 (fx1+ i)]
        [l (list) (if (vector-ref vec i) (cons i l) l)])
      ((fx=? i (vector-length vec)) l)))
                                                                                            ;
(define (vector-refs vec refs)           ;; (Vectorof X) (Vectorof Nat) -> (Vectorof X)
  ;; produce vector of selected elements of vec indexed by refs
  (let ([vrefs (make-vector (vector-length refs))])
    (do ([i 0 (fx1+ i)]) ((fx=? i (vector-length refs)) vrefs)
      (vector-set! vrefs i (vector-ref vec (vector-ref refs i))))))
                                                                                            ;
(define (list->bitwise ns)               ;; (Listof Nat) -> Bits
  ;; produce bitwise value by setting bits indexed by elements of ns
  (fold-left (lambda (acc n) (bitwise-copy-bit acc n 1)) 0 ns))
                                                                                            ;
(define (vector->bitwise vec)            ;; (Vectorof Nat) -> Bits
  (list->bitwise (vector->list vec)))
                                                                                            ;
(define (bitwise-span bits)              ;; Bits -> Nat
  ;; produce span of 1 bits in bits
  (fx- (bitwise-length bits) (bitwise-first-bit-set bits)))
                                                                                            ;
(define (bitwise->list bits)             ;; Bits -> (Listof Nat)
  ;; produce indices in ascending order of set bits in bits
  ;; negative bits (infinite 1s) is clipped to fixnum range!
  (let loop (
        [bits (if (negative? bits)
                (bitwise-and bits (greatest-fixnum))
                bits)]
        [result (list)])
    (if (zero? bits)  (reverse! result)
        (let ([b (bitwise-first-bit-set bits)])
          (loop (bitwise-copy-bit bits b 0) (cons b result))))))
                                                                                            ;
(define (vector-sample source size)      ;; (Vectorof X) Nat -> (Vectorof X)
  ;; produce random selection of length size from source (using Durstenfeld shuffle)
  (let* ( [source-length (vector-length source)]
          [size    (fxmin size source-length)]
          [shuffle (if (fx<? size source-length) size (fx1- size))])
    (cond [(or (fxzero? size) (fxzero? source-length))  '#()]
          [(fx=? size 1)
            (vector (vector-ref source (random source-length)))]
          [else
            (let ([vec (vector-take source-length source)])
              (do ([n 0 (fx1+ n)]) ((fx=? n shuffle) (vector-take size vec))
                (let* ( [r (fx+ n (random (fx- source-length n)))]
                        [t (vector-ref vec r)])
                  (vector-set! vec r (vector-ref vec n))
                  (vector-set! vec n t))))])))
                                                                                            ;
(define (u32-sample source size)         ;; {U32} Nat -> {U32}
  ;; produce random selection of length size from source (using Durstenfeld shuffle)
  (cond
    [(or (fxzero? size) (null? source))  '()]
    [(fx=? size 1)
       (list (list-ref source (random (length source))))]
    [else
      (let* ( [source (uint-list->bytevector source (native-endianness) 4)]
              [source-length (bytevector-length source)]
              [size   (fxmin (fx* size 4) source-length)])
        (do ([n 0 (fx+ n 4)])
            ((fx=? n size)
              (bytevector-truncate! source size)
              (bytevector->uint-list source (native-endianness) 4))
          (let* ([r (fx+ n (fxand #xFFFFFFFC (random (fx- source-length n))))]
                 [t (bytevector-u32-native-ref source r)])
            (bytevector-u32-native-set! source r (bytevector-u32-native-ref source n))
            (bytevector-u32-native-set! source n t)))) ] ))
                                                                                            ;
(define (fxsearch v target)              ;; (FXVectorof Fixnum) Fixnum -> Fixnum | #f
  ;; binary search [no benefit from using bytevectors or unroll?]
  (let search ([left 0] [right (fx1- (fxvector-length v))])
    (and (fx<=? left right)
      (let* ( [mid     (fxdiv (fx+ left right) 2)]
              [element (fxvector-ref v mid)]) 
        (cond 
          [ (fx<? element target) (search (fx1+ mid) right     ) ]
          [ (fx<? target element) (search left       (fx1- mid)) ]
          [else element])))))
                                                                                            ;
#;(define-syntax fxsearch                ;; (FXVectorof Fixnum) Fixnum -> Fixnum | #f
  ;; binary search (force inlined)
  (lambda (x) (syntax-case x ()
  [ (fxsearch v target)
    #'(let search ([left 0] [right (fx1- (fxvector-length v))])
        (and (fx<=? left right)
          (let* ( [mid     (fxdiv (fx+ left right) 2)]
                  [element (fxvector-ref v mid)]) 
            (cond 
              [ (fx<? element target) (search (fx1+ mid) right) ]
              [ (fx<? target element) (search left  (fx1- mid)) ]
              [else element]))))])))
                                                                                            ;
(define (key-word-args args defaults)    ;; {KWarg} {KWarg} -> {X}
  ;; produce values from defaults overridden by args values with matching key
  (map (lambda (default-kv)
      (let ([kv (assq (car default-kv) args)])  ;; if this default in args
        (if kv (cdr kv) (cdr default-kv))))     ;; then use given val
    defaults))
                                                                                            ;
(define (do-with-progress n f)           ;; Nat (Nat -> ) ->
  ;; apply f to 0..n-1 with display of iteration time
  (define (secs-since t) 
    (quotient (+ (- (cpu-time) t) 500) 1000))
  (display " starting     \r")
  (flush-output-port (current-output-port))
  (let* ( [start-t0 (cpu-time)]
          [stride   #;(expt 10 (exact (ceiling (log (max 1 (log (fxmax 1 n)))))))
                    1 #;(isqrt (fxmax 1 n))])
    (do ([step 0 (fx+ step stride)]) ((fx>=? step n))
      (let* ([limit   (fxmin n (fx+ step stride))]
             [step-t0 (cpu-time)])
        (do ([i step (fx1+ i)]) ((fx=? i limit))
          (f i))
        (for-each display
          `(#\space ,limit ": " ,(secs-since step-t0) "  " ,(secs-since start-t0) "     \r"))
        (flush-output-port (current-output-port))))))
                                                                                            ;
(define (intersect1d l1 l2)              ;; {Fixnum} {Fixnum} -> {Fixnum}
  ;; produce intersection of sorted lists of fixnums, without duplicates
  ;  (assert (equal? l1 (list-sort fx<? l1)))
  ;  (assert (equal? l2 (list-sort fx<? l2)))
  ;  (let ((result
  (let loop ([l1 l1] [l2 l2] [result (list)])
    (if (pair? l1)
      (let loop2 ([carl1 (car l1)] [l2 l2] [result result])
        (if (pair? l2)
          (cond [(fx<? carl1 (car l2)) 
                    (loop (cdr l1) l2 result)]
                [(fx>? carl1 (car l2)) 
                    (loop2 carl1 (cdr l2) result)]
                [else (loop (cdr l1) (cdr l2) 
                        (if (and (pair? result) (fx=? carl1 (car result)))
                          result
                          (cons carl1 result)))])
          (reverse! result)))
      (reverse! result))))
  ;  ) (assert (equal? result (unique! fx=? (append result '()))))
  ;    (assert (equal? result (list-sort fx<? result)))
  ;    result))
                                                                                            ;
(define union1d                          ;; {Fixnum} {Fixnum} [ {Fixnum} ] -> {Fixnum}
  ;; produce union of sorted lists of fixnums [(merge fx<? l1 l2) without duplicates]
  (case-lambda
    [ (ls1 ls2)
      (unique! fx=? (cond
        [(null? ls1) ls2]
        [(null? ls2) ls1]
        [(fx<? (car ls2) (car ls1))
          (cons (car ls2) (union1d ls1 (cdr ls2)))]
        [else (cons (car ls1) (union1d (cdr ls1) ls2))])) ]
    [ (ls1 ls2 ls3)
      (unique! fx=? (cond
        [(null? ls1) (union1d ls2 ls3) ]
        [(null? ls2) (union1d ls1 ls3) ]
        [(null? ls3) (union1d ls1 ls2) ]
        [else
          (let ([min (fxmin (car ls1) (car ls2) (car ls3))])
            (cons min (union1d (if (fx=? min (car ls1)) (cdr ls1) ls1)
                               (if (fx=? min (car ls2)) (cdr ls2) ls2)
                               (if (fx=? min (car ls3)) (cdr ls3) ls3)))) ])) ]))
                                                                                            ;
(define (setdiff1d l1 l2)                ;; {Fixnum} {Fixnum} -> {Fixnum}
  ;; produce difference of sorted lists of fixnums
  ;  (assert (equal? l1 (list-sort fx<? l1)))
  ;  (assert (equal? l2 (list-sort fx<? l2)))
  ;  (let ((result
  (let loop ([l1 l1] [l2 l2] [result (list)])
    (cond [(null? l1)  (reverse! result)]
          [(or (null? l2) (fx<? (car l1) (car l2)))
              (loop (cdr l1) l2 (cons (car l1) result))]
          [(fx>? (car l1) (car l2))
              (loop l1 (cdr l2) result)]
          [else (loop (cdr l1) (cdr l2) result)])))
  ;  ) (assert (equal? result (unique! fx=? (append result '()))))
  ;    (assert (equal? result (list-sort fx<? result)))
  ;    result))
                                                                                              ;
(define (in1d x1s x2s)                   ;; {X} {X} -> Bits
  ;; produce index mask of x1s that are also in x2s
  (let next-x1 ([x1s x1s] [bits 0] [bitx 0])
    (if (null? x1s)  bits
        (next-x1 (cdr x1s)
                 (if (memv (car x1s) x2s)
                   (bitwise-copy-bit bits bitx 1)
                   bits)
                 (fx1+ bitx)))))
                                                                                            ;
(define (include-by-mask xs mask)        ;; {X} Bits -> {X}
  ;; extract elements of xs corresponding to 1 bits in mask
  (vector->list
    (vector-refs 
      (list->vector xs) 
      (list->vector (bitwise->list mask)))))
                                                                                            ;
(define (exclude-by-mask xs mask)        ;; {X} Bits -> {X}
  ;; omit elements of xs corresponding to 1 bits in mask
  (cond
    [(zero? mask) xs]
    [(odd?  mask)
      (exclude-by-mask (cdr xs) (bitwise-arithmetic-shift-right mask 1))]
    [else
      (cons (car xs) (exclude-by-mask (cdr xs) (bitwise-arithmetic-shift-right mask 1)))]))
                                                                                            ;
  (define thread-limit 7)                ;; #hyperthreads-1 for best wall time, #cores for best cpu?
                                                                                            ;                                                                                            ;
  (define thread-if 3)                   ;; fork if vector-length >= thread-if
                                                                                            ;                                                                                            ;
  (define random-seed-limit (- (expt 2 32) 1))
                                                                                            ;                                                                                            ;
  (define seed (make-thread-parameter #f))
                                                                                            ;                                                                                            ;
(define (threaded-vector-for-each f . vs);; Nat (X ... -> ) (Vectorof X) ... ->
  ;; in a new thread for each, apply f to elements of vs, return when all finished
  (let ([todo (vector-length (car vs))])
    (if (fx<? todo thread-if)
      (apply vector-for-each (lambda xs  ;; could be (apply vector-for-each f vs)]
          (apply f xs))
        vs)
      (let ([threads  0]
            [mutex   (make-mutex)]
            [free    (make-condition)]
            [done    (make-condition)])
        (apply vector-for-each (lambda xs
            (with-mutex mutex
              (when (fx>=? threads thread-limit)
                (condition-wait free mutex))
              (set! threads (fx1+ threads)))
            (seed (fx1+ (random random-seed-limit)))
            (fork-thread (lambda ()
                    (random-seed (seed))               ;; child thread
                    (apply f xs)                       ;;
                    (with-mutex mutex                  ;;
                      (set! todo (fx1- todo))          ;;
                      (when (fxzero? todo)             ;;
                        (condition-signal done))       ;;
                      (set! threads (fx1- threads))    ;;
                      (condition-signal free)))        ;;
              ))
          vs)
        (with-mutex mutex
          (unless (fxzero? todo)
            (condition-wait done mutex)))))))
                                                                                            ;
(define cost-center-1 (make-cost-center))
(define cost-center-2 (make-cost-center))
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
  [expect ( [unique! fx=? '(1 2 2 3 4 4)] '(1 2 3 4) )]
  [expect ( [condense! '(1 2 1 1 3 2 3 1)] '(1 2 3) )]
  [expect ( [argmax-multi '(2 1 3 4 -1) '(2 2 1 1 3)]  '(0 3 4))]                                                                                           
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
  [expect ( [vector-fold-left (lambda (l x) (cons x l)) '() '#(0 1 2)] '(2 1 0) )
          ( [vector-fold-left (lambda (s x y) (* s (+ x y))) 1 '#(0 1 2) '#(1 2 3)] 15 )]
  [expect ( [vector-average '#(1 2 5)] 8/3 )]
  [expect ( [fxvector-max '#(0 -3 7 5)] 7 )]
  [expect ( [fx3<- '#(1/3 1.0)] '#(333 1000) )]
  [expect ( [fx3* 500 500] 250 )]
  [expect ( [vector-take 0 '#(0 1 2 3)] '#()       )
          ( [vector-take 2 '#(0 1 2 3)] '#(0 1)    )
          ( [vector-take 5 '#(0 1 2 3)] '#(0 1 2 3))]
  [expect ( [vector-count zero? '#(0 1 2 0 3)] 2 )
          ( [vector-count not   '#(#t 0 fx1+)] 0 )]
  [expect ( [vector-indices '#(#f #f #f)] '()    )
          ( [vector-indices '#(#t #f  0)] '(2 0) )]
  [expect ( [vector-refs '#(0 11 22 33 44) '#(0 2 4)] '#(0 22 44) )]
  [expect ( [bitwise-span #b01011101000] 7 )]
  (let* ( [source (build-vector 100 id)]
          [select (vector->list [vector-sample source 50])])
    [expect ( [length select] 50 )
            ( [list? (for-all memq select (make-list 50 (vector->list source)))] #t )  ;; all from source
            ( [map (lambda (x) (length (remq x select))) select] (make-list 50 49)) ]) ;; all different
  [expect ( [(lambda (x . args)          ;; example fn with kws k1 & k2: default is (x 1 2)
              (append (list x) [key-word-args args '([k1 . 1] [k2 . 2])] ))
            99 '[k2 . 22] '[k1 . 11] '[k3 . 33] ] '(99 11 22) )]           ;; apply fn to 99 '[k2 ...
  [expect ( [intersect1d     '(1 2 3 4) '(1 3 5)] '(1 3) )]
  [expect ( [intersect1d     '(1 1 2)   '(1 1 3)] '(1) )]
  [expect ( [union1d         '(1 2 4 5) '(1 3 5)] '(1 2 3 4 5) )]
  [expect ( [union1d         '(1 1 2)   '(1 1 3)] '(1 2 3) )]
  [expect ( [union1d         '(1 4) '(2 5 5) '(3 6)] '(1 2 3 4 5 6) )]
  [expect ( [setdiff1d       '(1 2 3 4) '(1 3 5)] '(2 4) )]
  [expect ( [in1d            '(1 2 3 4) '(1 3 5)]  #b101 )]
  [expect ( [include-by-mask '(1 2 3 4)  #b1101]  '(1 3 4))]
  [expect ( [exclude-by-mask '(1 2 3 4)  #b1101]  '(2))]
  
  #;(let ([xs (make-vector 3 #f)])
    [expect ( [begin [threaded-vector-for-each (lambda (x ix)
                          (vector-set! xs ix (+ x 1)))
                        '#(1 2 3) '#(0 1 2)]
                      xs]   '#(2 3 4))])
  
)