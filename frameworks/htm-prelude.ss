;; © 2019 Roger Turner <https://github.com/rogerturner/HTM-scheme/issues/new/choose>
;; SPDX-License-Identifier: AGPL-3.0-or-later  (see Notices below)

#| HTM-scheme Prelude

Common functions specialized for HTM-scheme (Fixnums used wherever possible,
includes mutating and proper-list assuming versions of standard procedures).

See htm-concept.ss for an overviw of system objectives and code conventions.
Indentation facilitates using an editor "Fold All" view for a file overview.

Functions are documented by type and purpose comments, and checked examples.

Types:
  Boolean, Number, Integer, Fixnum, (Listof X), (Vectorof X) ... = Scheme types
  X Y -> Z = function of arg types X, Y to result type Z, []=optional, !=mutated
  {X}      = abbreviation for (Listof X)
  T?       = truthy (everything except #f)
  Nat      = natural number (non-negative Fixnum or exact Integer)
  Fixnum3  = Fixnum interpreted as number with 3 decimal places, ie 1000 = 1.0
  Bits     = Integer (Bignum) interpreted bitwise
  KWarg    = Pair (key . value)
  
"Lisp's parentheses are the bumps on the top of Lego" [Paul Graham]
"car and cdr are the only honest function names"  [citation needed]

  |#

  #!chezscheme

(library (frameworks htm-prelude)
                                                                                            ;
(export
  *examples*
  example:
  check-examples
  build-list
  take
  list-average
  !reverse!
  !partition
  list-refs!
  unique!
  except-last!
  flattenr
  filter!
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
  fxvector-fold-left
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
  list->u32-vector
  u32-search
  u32-sample
  list-sample
  fxmemv
  fxsearch
  string->word-list
  bytevector-bit-set?
  bytevector-copy-bit!
  bytevector-copy-bit
  key-word-args
  display-args
  do-with-progress
  intersect1d
  union1d
  union1d!
  setdiff1d
  setdiff1d!
  in1d
  include-by-mask
  exclude-by-mask
  make-thread-pool
  threaded-for-each-x
  threaded-vector-for-each
  cost-center-1
  with-thunk-time
  reset-thunk-time
  display-thunk-time
  vr0@
  vr0!
  vr1@
  vr1!
  first
  define-memoized
  expect
  expect-count
  expects
  )
                                                                                            ;
(import (chezscheme))
                                                                                            ;
  (implicit-exports #f)

(meta define *examples* (list))          ;; {(SyntaxObject . SyntaxObject)}
                                                                                            ;
(meta define *report-examples* #f)       ;; Boolean
                                                                                            ;
(define-syntax example:                  ;; Expr | Expr => Result ->
  ;; (example: e) or (example: e => r)
  ;; cons syntax objects for e, 'e to *examples*; produce dummy syntax object;
  ;; e and r will be evaluated by check-examples, so cannot reference lexical environment;
  ;; crude implementation, but short and allows readable inline checked examples
  (lambda (e-expr)
    (set! *examples*
      (cons
        (syntax-case e-expr (=>)
          [(example: e)      (cons #'e #''e) ]
          [(example: e => r) (cons #'(equal? e r) #''e) ])
        *examples*))
    #'(define _ #f)))
                                                                                            ;
(define-syntax check-examples            ;; -> (begin (assert check) ...)
  ;; produce asserts from *examples*, clear *examples*, maybe display number checked
  (lambda (_)
    (let ([examples (reverse *examples*)])
      (set! *examples* (list))
      (syntax-case _ ()
        [(_) #`(begin
              #,@(map (lambda (ex)
                   (let ([c (car ex)] [e (cdr ex)])
                    #`(unless #,c
                        (assertion-violation #f "example" #,e))))
                  examples)
                (when #,*report-examples*
                  (display #,(length examples))
                  (display " examples checked\n"))) ]))))
  
  (example: (+ 2 2) => 4 )
                                                                                            ;
  (define (double n)                     ;; example in definition
    (example: (double 2) => 4 )
    (+ n n))
    
  (example:                              ;; example with custom test condition
    (let-values ([(q r) (div-and-mod 17.5 3)])
      (and (eqv? q 5.0) (eqv? r 2.5))) )
                                                                                            ;
(define (build-list n f)                 ;; Nat (Nat -> X) -> (Listof X)
  ;; produce list of n X's by applying f to 0..n-1 (in reverse order)
  (example: (build-list 3 -) => '(0 -1 -2) )
  (let next-x ([i (fx1- n)] [xs (list)])
    (if (fxnegative? i)  xs
        (next-x (fx1- i) (cons (f i) xs)))))
                                                                                            ;
(define (take n xs)                      ;; Nat (Listof X) -> (Listof X)
  ;; produce first n elements of xs
  (example: (take 2 '(a 2 c)) => '(a 2) )
  (example: (take 4 '(a 2 c)) => '(a 2 c) )
  (list-head xs (fxmin n (length xs))))
                                                                                            ;
(define (list-average ns)                ;; (Listof Number) -> Number
  ;; produce mean of xs (a non-empty list)
  (example: (list-average '(1 2 27)) => 10 )
  (/ (apply + ns) (length ns)))
                                                                                            ;
(define (!reverse! xs)                   ;; {X}! -> {X}
  ;; mutate xs (a proper list) to reverse
  (example: (!reverse! (list 1 2 3)) => '(3 2 1))
  (example: (let* ([xs (list 1 2)] [x2 (cdr xs)]) (eq? (!reverse! xs) x2)))
  (let loop ([xs xs] [rs (list)])
    (cond
      [(null? xs) rs ]
      [else (let ([next (cdr xs)])
              (set-cdr! xs rs)
              (loop next xs)) ])))
                                                                                            ;
(define (!partition pred? xs)            ;; (X -> T?) {X} -> {X} {X}
  ;; produce partition of xs (a proper list) by pred?
  (example: (let-values ([(ts fs) (!partition odd? '(1 2 3))])
              (list ts fs)) => '((1 3) (2)) )
  (let loop ([xs xs] [ts (list)] [fs (list)])
    (cond
      [(null? xs) (values (!reverse! ts) (!reverse! fs)) ]
      [(pred? (car xs)) (loop (cdr xs) (cons (car xs) ts) fs) ]
      [else (loop (cdr xs) ts (cons (car xs) fs)) ])))
                                                                                            ;
(define (list-refs! ls xs)               ;; {X}! {Nat} -> {X}
  ;; mutate ls by removing elements not in xs, which are sorted indexes of ls elements
  (example: (list-refs! (list 0 1 2 3) '(1 3)) => '(1 3))
  (if (null? xs)  '()
    (let ([ls (list-tail ls (car xs))])
      (let next-x ([prev ls] [xs (cdr xs)] [prev-x (car xs)])
        (if (null? xs)  (begin (set-cdr! prev '())  ls)
          (let ([this-x (car xs)])
            (unless (fx=? this-x (fx1+ prev-x))
              (set-cdr! prev (list-tail prev (fx- this-x prev-x))))
            (next-x (cdr prev) (cdr xs) this-x)))))))
                                                                                            ;
(define (unique! eql? xs)                ;; (X X -> T?) {X}! -> {X}
  ;; mutate xs removing adjacent duplicates by eql?
  (example: (unique! = (list 1 2 2 3 4 4)) => '(1 2 3 4) )
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
(define (except-last! xs)                ;; {X}! -> {X}
  ;; produce xs without last element
  (example: (except-last! '(1 2 3)) => '(1 2) )
  (if (or (null? xs) (null? (cdr xs)))  '()
      (let next-x ([p xs])
        (cond
          [(null? (cddr p))
            (set-cdr! p '())
            xs ]
          [else (next-x (cdr p)) ]))))
                                                                                            ;
(define (flattenr xss)                   ;; { {X} } -> {X}
  ;; produce list of all x (reversed)
  (example: (flattenr '((1 2) () (3))) => '(3 2 1) )
  (let next-list ([xss xss] [result (list)])
    (cond
      [(null? xss) result ]
      [else
        (let next-xs ([xs (car xss)] [result result])
          (cond
            [(null? xs) (next-list (cdr xss) result) ]
            [else (next-xs (cdr xs) (cons (car xs) result)) ])) ])))
                                                                                            ;
(define (filter! pred? xs)               ;; (X -> T?) {X}! -> {X}
  ;; produce mutated xs, omitting elements for which (pred? x) is false
  (example: (filter! odd? (list 1 2 3 4 4 5 5)) => '(1 3 5 5) )
  (example: (filter! odd? (list 0 0 1 3 5 6 6)) => '(1 3 5) )
  (example: (filter! odd? (list 0 2 4 6))       => '() )
  (example: (let ([xs (list 0 1)]) (eq? (filter! odd? xs) (cdr xs) )))
  (let next-head ([head xs])
    (cond
      [(null? head) head ]
      [(pred? (car head))
        (let next-in ([last head] [rest (cdr head)])
          (cond
            [(null? rest)
              (set-cdr! last rest)
              head ]
            [(pred? (car rest))
              (next-in rest (cdr rest)) ]
            [else (let next-out ([next (cdr rest)])
                (cond
                  [(null? next)
                    (set-cdr! last next)
                    head ]
                  [(pred? (car next))
                    (set-cdr! last next)
                    (next-in next (cdr next)) ]
                  [else (next-out (cdr next)) ])) ])) ]
      [else (next-head (cdr head)) ])))
                                                                                            ;
(define (condense! xs)                   ;; {X}! -> {X}
  ;; mutate xs by removing eq? duplicates
  (example: (condense! (list 1 2 1 1 3 2 3 1)) => '(1 2 3) )
  (let condense ([len (length xs)] [xs xs])
    (if (fx<? len 2)  xs
      (let* ( [half-len (fxdiv len 2)]
              [last (list-tail xs (fx1- half-len))]
              [rest (cdr last)])
        (set-cdr! last '())
        (let* ( [filtered-begin (condense half-len xs)]
                [filtered-end
                  (let ([xs (filter! (lambda (x)
                               (not (memq x filtered-begin)))
                             rest)])
                    (condense (length xs) xs))])
          (append! filtered-begin filtered-end))))))
                                                                                            ;
(define (argmax-multi a group-keys)      ;; {Fixnum} {Fixnum} -> {Nat}
  ;; Get the indices of the max values of each group in 'a', grouping the
  ;; elements by their corresponding value in group-keys
  (example: (argmax-multi '(2 1 3 4 -1) '(2 2 1 1 3)) => '(0 3 4) )                                                                                           
  (let per-group ([a a] [group-keys group-keys] [index 0] [result (list)])
    (cond
      [ (null? group-keys) (!reverse! result) ]
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
  ;; produce the nearest integer to x
  (example: (int<- 2.500000000000001) => 3 )
  (example: (int<- 5/2              ) => 2 )
  (example: (int<- 7/2              ) => 4 )
  (exact (round x)))
                                                                                            ;
(define (id x)                           ;; X -> X
  x)
                                                                                            ;
(define (build-vector n f)               ;; Nat (Nat -> X) -> (Vectorof X)
  ;; produce vector length n by applying f to indices
  (example: (build-vector 3 id) => '#(0 1 2) )
  (let ([v (make-vector n)])
    (do ([i 0 (fx1+ i)]) ((fx=? i n) v)
      (vector-set! v i (f i)))))
                                                                                            ;
(define (vector-filter pred vec)         ;; (X -> T?) (Vectorof X) -> (Vectorof X)
  ;; produce vector of elements of vec for which (pred elt) is not false
  (example: (vector-filter even? '#(1 2 3 4 5)) => '#(2 4) )
  (list->vector (filter! pred (vector->list vec))))
                                                                                            ;
(define (indexes seq)                    ;; (Vectorof X) -> (Vectorof Nat)                                                                 
  ;; produce indexes of vector or list   ;; (Listof X)   -> (Listof Nat)
  (example: (indexes '(a b c))  => '(0 1 2)  )
  (example: (indexes '#(1 2 3)) => '#(0 1 2) )
  (if (vector? seq)
    (build-vector (vector-length seq) id)
    (iota         (length seq))))
                                                                                            ;
(define (with-index folder f vs)         ;; ((? -> ?) ?... -> ?) (?... -> ?) (Listof (Vectorof ?)) -> ?
  (apply folder
    (lambda xs (apply f xs))
    (append vs (list (indexes (car vs))))))
                                                                                            ;
(define (vector-map-x f . vs)            ;; (X ... Nat -> Y) (Vectorof X) ... -> (Vectorof Y)
  ;; produce vector by applying f to each element of vs and its index
  (example: (vector-map-x (lambda (x y i) (+ x y i)) '#(10 11 12) '#(1 1 1)) => '#(11 13 15) )
  (with-index vector-map f vs))
                                                                                            ;
(define (vector-for-each-x f . vs)       ;; (X ... Nat -> ) (Vectorof X) ... -> 
  ;; apply f to each element of vs and its index
  (with-index vector-for-each f vs))
                                                                                            ;
(define (vector-fold-left f o . vs)      ;; (X Y ... -> X) X (Vectorof Y) ... -> X
  ;; assume all vs same length
  (example: (vector-fold-left (lambda (l x) (cons x l)) '() '#(0 1 2)) => '(2 1 0) )
  (example: (vector-fold-left (lambda (s x y) (* s (+ x y))) 1 '#(0 1 2) '#(1 2 3)) => 15 )
  (let ([lastx (fx1- (vector-length (car vs)))])
    (if (fxnegative? lastx)  o
        (let fold ([acc o] [vx 0])
          (define (do-apply)
            (let vs/x->list ([vs vs] [xs (list)])
              (if (null? vs)  (apply f acc (reverse xs))
                  (vs/x->list (cdr vs) (cons (vector-ref (car vs) vx) xs)))))
          (if (fx=? vx lastx)  (do-apply)
            (fold (do-apply) (fx1+ vx)))))))
                                                                                            ;
(define (fxvector-fold-left f o v)       ;; (X Y -> X) X (FXVectorof Y) -> X
  (example: (fxvector-fold-left (lambda (s x) (+ s (* x x))) 0 '#vfx(1 2 3)) => 14 )
  (let ([lastx (fx1- (fxvector-length v))])
    (if (fxnegative? lastx)  o
        (let fold ([acc o] [vx 0])
          (if (fx=? vx lastx)
            (f acc (fxvector-ref v vx))
            (fold (f acc (fxvector-ref v vx)) (fx1+ vx)))))))
                                                                                            ;
(define (vector-average vec)             ;; (Vectorof Number) -> Number
  ;; (list-average (vector->list vec))
  (example: (vector-average '#(1 2 5)) => 8/3 )
  (/ (vector-fold-left + 0 vec) (vector-length vec)))
                                                                                            ;
(define (fxvector-max vec)               ;; FXVector -> Fixnum
  ;; produce max element of vec
  (example: (fxvector-max '#vfx(0 -3 7 5)) => 7 )
  (fxvector-fold-left fxmax (least-fixnum) vec))
                                                                                            ;
(define fx3 1000)                        ;; Fixnum                                                                                            
                                                                                            ;
(define (fx3<- x)                        ;; [Vectorof] Number -> [Vectorof] Fixnum3
  ;; produce Fixnum3 representation of number or vector of numbers
  (example: (fx3<- '#(1/3 1.0)) => '#(333 1000) )
  (if (vector? x) (vector-map fx3<- x)
                  (int<- (* x fx3))))
                                                                                            ;
(define (fx3* x y)                       ;; Fixnum3 Fixnum3 -> Fixnum3
  ;; produce product
  (example: (fx3* 500 500) => 250 )
  (fxdiv (fx+ (fx* x y) (fxdiv fx3 2)) fx3))
                                                                                            ;
(define (vector-extend vec . incr)       ;; (Vectorof X) [Nat] -> (Vectorof X)
  ;; increase the length of vec by incr, or double it if no incr
  (let* ( [len    (vector-length vec)]
          [result (make-vector
                    (if (null? incr)  (fx* 2 len)
                        (fx+ len (car incr))))])
    (do ([i 0 (fx1+ i)]) ((fx=? i len) result)
      (vector-set! result i (vector-ref vec i)))))
                                                                                            ;
(define (vector-take n vec)              ;; Nat (Vectorof X) -> (Vectorof X)
  ;; produce copy of first n elements of vec, or copy of vec if n > length
  (example: (vector-take 0 '#(0 1 2 3)) => '#()       )
  (example: (vector-take 2 '#(0 1 2 3)) => '#(0 1)    )
  (example: (vector-take 5 '#(0 1 2 3)) => '#(0 1 2 3))
  (let* ( [size   (fxmin n (vector-length vec))]
          [result (make-vector size)])
    (do ([i 0 (fx1+ i)]) ((fx=? i size) result)
      (vector-set! result i (vector-ref vec i)))))
                                                                                            ;
(define (vector-count proc vec)          ;; (X -> T?) (Vectorof X) -> Nat
  ;; produce count of elements of vec for which (proc elt) is not false
  (example: (vector-count zero? '#(0 1 2 0 3)) => 2 )
  (example: (vector-count not   '#(#t 0 fx1+)) => 0 )
  (vector-fold-left (lambda (acc elt)
      (if (proc elt)  (fx1+ acc)  acc ))
    0 vec))
                                                                                            ;
(define (vector-indices vec)             ;; (Vectorof X) -> (Listof Nat)
  ;; produce list of indices for which vector element is not false
  (example: (vector-indices '#(#f #f #f)) => '()    )
  (example: (vector-indices '#(#t #f  0)) => '(0 2) )
  (do ( [i 0 (fx1+ i)]
        [l (list) (if (vector-ref vec i)  (cons i l)  l)])
        ((fx=? i (vector-length vec)) (!reverse! l))))
                                                                                            ;
(define (vector-refs vec refs)           ;; (Vectorof X) (Vectorof Nat) -> (Vectorof X)
  ;; produce vector of selected elements of vec indexed by refs
  (example: (vector-refs '#(0 11 22 33 44) '#(0 2 4)) => '#(0 22 44) )
  (let* ( [nrefs (vector-length refs)]
          [vrefs (make-vector nrefs)])
    (do ([i 0 (fx1+ i)]) ((fx=? i nrefs) vrefs)
      (vector-set! vrefs i (vector-ref vec (vector-ref refs i))))))
                                                                                            ;
(define (list->bitwise ns)               ;; (Listof Nat) -> Bits
  ;; produce bitwise value by setting bits indexed by elements of ns
  (example: (list->bitwise '(0 2 3)) => #b1101 )
  (fold-left (lambda (acc n) (bitwise-copy-bit acc n 1)) 0 ns))
                                                                                            ;
(define (vector->bitwise vec)            ;; (Vectorof Nat) -> Bits
  (list->bitwise (vector->list vec)))
                                                                                            ;
(define (bitwise-span bits)              ;; Bits -> Nat
  ;; produce span of 1 bits in bits
  (example: (bitwise-span #b01011101000) => 7 )
  (fx- (bitwise-length bits) (bitwise-first-bit-set bits)))
                                                                                            ;
(define (bitwise->list bits)             ;; Bits -> (Listof Nat)
  ;; produce indices in ascending order of set bits in bits
  ;; negative bits (infinite 1s) is clipped to fixnum range!
  (example: (bitwise->list #b1101) => '(0 2 3))
  (let loop (
        [bits (if (negative? bits)
                (bitwise-and bits (greatest-fixnum))
                bits)]
        [result (list)])
    (if (zero? bits)  (!reverse! result)
        (let ([b (bitwise-first-bit-set bits)])
          (loop (bitwise-copy-bit bits b 0) (cons b result))))))
                                                                                            ;
(define (vector-sample source size)      ;; (Vectorof X) Nat -> (Vectorof X)
  ;; produce random selection of length size from source (using Durstenfeld shuffle)
  (example: (let ([s (vector-sample '#(0 1) 1)])
              (or (equal? s '#(0)) (equal? s '#(1)))))
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
  (alias fxasl  fxarithmetic-shift-left)
  (alias fxasr  fxarithmetic-shift-right)
  (alias bvu32@ bytevector-u32-native-ref)
  (alias bvu32! bytevector-u32-native-set!)
                                                                                            ;
(define (list->u32-vector u32s)          ;; {Fixnum} -> Bytevector
  ;; produce low 32 bits of u32s
  (example: (list->u32-vector '(1 -2)) => '#vu8(1 0 0 0 254 255 255 255) )
  (let* ( [n  (fx* 4 (length u32s))]
          [bv (make-bytevector n)])
    (do ( [i 0 (fx+ i 4)]
          [u32s u32s (cdr u32s)])
        ((fx=? i n) bv)
      (bvu32! bv i (fxand #xFFFFFFFF (car u32s))))))
                                                                                            ;
(define (u32-search u32s target)         ;; Bytevector Fixnum -> Fixnum | #f
  ;; produce target if in u32s or #f by binary search
  (example: (u32-search '#vu8(1 0 0 0 2 0 0 0 3 0 0 0) 3) => 3 )
  (example: (u32-search '#vu8(1 0 0 0 2 0 0 0 3 0 0 0) 0) => #f )
  (let search ([left 0] [right (fx- (bytevector-length u32s) 4)])
    (and (fx<=? left right)
      (let* ( [mid  (fxasl (fxasr (fx+ left right) 3) 2)]
              [u32  (bvu32@ u32s mid)])
        (cond 
          [ (fx<? u32 target) (search (fx+ mid 4) right)]
          [ (fx<? target u32) (search left (fx- mid 4)) ]
          [ else u32 ])))))
                                                                                            ;
(define (u32-sample source size)         ;; {U32} Nat -> {U32}
  ;; produce random selection of length size from source (using Durstenfeld shuffle)
  (example: (let ([s (u32-sample '(0 1) 1)])
              (or (equal? s '(0)) (equal? s '(1)))))
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
                 [t (bvu32@ source r)])
            (bvu32! source r (bvu32@ source n))
            (bvu32! source n t)))) ] ))
                                                                                            ;
  (example: (let ([s (list-sample '(0 1) 1)])
              (or (equal? s '(0)) (equal? s '(1)))))
(define list-sample                      ;; {X} Nat [Nat] -> {X}
  ;; produce random selection (in order) of length size from xs (length n-xs > 0)
  ;; (not unbiased, but close enough with HTM cellx lists?)
  (case-lambda
    [(xs size)
      (list-sample xs size (length xs)) ]
    [(xs size n-xs)
      (cond
        [(fxnonpositive? size)  '() ]
        [(fx=? size 1)  (list (list-ref xs (random n-xs))) ]
        [else
          (let ([stride (fxmin (fx- n-xs size) (fxdiv (fx* 2 n-xs) size))])
            (let sample ([xs xs] [n size] [n-xs n-xs] [skip (random stride)] [out (list)])
              (cond
                [(fxzero? skip)
                  (let ([out (cons (car xs) out)]
                        [n   (fx1- n)])
                    (if (fxzero? n)  (!reverse! out)
                      (let* ( [n-xs   (fx1- n-xs)]
                              [stride (fx1+ (fxdiv n-xs size))])
                        (sample (cdr xs) n n-xs (random stride) out)))) ]
                [else  
                  (sample (cdr xs) n (fx1- n-xs) (fx1- skip) out) ]))) ]) ]))
                                                                                            ;
(define (fxmemv x xs)                    ;; Fixnum {Fixnum} -> {Fixnum} | #f
  ;; produce first tail of ls with car equal to x, or #f
  (example: (fxmemv 2 '(1 2 3)) => '(2 3) )
  (example: (fxmemv 4 '(1 2 3)) => #f     )
  (let fxmemv ([xs xs])
    (and (not (null? xs))
         (if (fx=? (car xs) x)  xs
             (fxmemv (cdr xs))))))
                                                                                            ;
(define (fxsearch v target)              ;; FXVector Fixnum -> Fixnum | #f
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
  #;
  (define-syntax fxsearch                ;; FXVector Fixnum -> Fixnum | #f
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
(define (string->word-list s)            ;; String -> {String}
  ;; produce list of words (delimited by space(s)) from s
  (example: (string->word-list "a  -b") => '("a" "-b") )
  (example: (string->word-list "a -b ") => '("a" "-b") )
  (example: (string->word-list "  ")    => '() ) 
  (let next-char ([s (string->list s)] [word (list)] [result (list)])
    (define (maybe-word)  ;; cons word to result
      (if (null? word) result
          (cons (list->string (reverse word)) result)))
    (cond
      [(null? s) (reverse (maybe-word)) ]
      [(char=? (car s) #\ )
        (next-char (cdr s) (list) (maybe-word)) ]
      [else (next-char (cdr s) (cons (car s) word) result) ])))
                                                                                            ;
(define (bytevector-bit-set? bv n)       ;; Bytevector Nat -> Boolean
  ;; produce #t if nth bit of bv is one, #f otherwise
  (example: (bytevector-bit-set? '#vu8(4) 2) => #t )
  (let-values ([(byte bit) (fxdiv-and-mod n 8)])
    (fxbit-set? (bytevector-u8-ref bv byte) bit)))
                                                                                            ;
(define (bytevector-copy-bit! bv n b)    ;; Bytevector! Nat (0 | 1) ->
  ;; set nth bit of bv to b
  (let-values ([(byte bit) (fxdiv-and-mod n 8)])
    (bytevector-u8-set! bv byte
      (fxcopy-bit (bytevector-u8-ref bv byte) bit b))))
                                                                                            ;
(define (bytevector-copy-bit bv n b)     ;; Bytevector Nat (0 | 1) -> Bytevector
  ;; produce copy of bv with nth bit set to b
  (let ([bv (bytevector-copy bv)])
    (bytevector-copy-bit! bv n b)
    bv))
                                                                                            ;
(define (key-word-args args defaults)    ;; {KWarg} {KWarg} -> {X}
  ;; produce values from defaults overridden by args values with matching key
  (map (lambda (default-kv)
      (let ([kv (assq (car default-kv) args)])  ;; if this default in args
        (if kv (cdr kv) (cdr default-kv))))     ;; then use given val
    defaults))
                                                                                            ;
(define (display-args args defaults)     ;; {KWarg} {KWarg} ->
  ;; display defaults overridden by args
  (for-each (lambda (default-kv)
      (for-each display `(
          "[" ,(car default-kv) " . "
          ,(let* ([kv (assq (car default-kv) args)]
                  [v  (if kv  (cdr kv) (cdr default-kv))])
              (cond
                [(pair? v) "(...)" ]
                [(vector? v) "#(...)" ]
                [else v]))
          "]\n")))
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
(define (intersect1d x1s x2s)            ;; {Fixnum} {Fixnum} -> {Fixnum}
  ;; produce intersection of sorted duplicate-free lists of fixnums
  (example: (intersect1d '(1 2 3 4) '(1 3 5)) => '(1 3) )
  (let loop ([x1s x1s] [x2s x2s] [out (list)])
    (if (or (null? x1s) (null? x2s))  (!reverse! out)
      (let ([x1 (car x1s)] [x2 (car x2s)])
        (cond
          [(fx<? x1 x2)  (loop (cdr x1s) x2s out) ]
          [(fx>? x1 x2)  (loop x1s (cdr x2s) out) ]
          [else  (loop (cdr x1s) (cdr x2s) (cons x1 out)) ])))))
                                                                                            ;
  (example: (union1d '(1 2 4 5) '(1 3 5)) => '(1 2 3 4 5) )
  (example: (union1d '(1 2 4) '(2 4 5) '(3 5 6)) => '(1 2 3 4 5 6) )
(define union1d                          ;; {Fixnum} {Fixnum} [ {Fixnum} ] -> {Fixnum}
  ;; produce union of sorted lists of fixnums; omit inter-list duplicates
  (case-lambda
    [ (ls1 ls2)
      (cond
        [(null? ls1) ls2]
        [(null? ls2) ls1]
        [else (let ([cls1 (car ls1)]
                    [cls2 (car ls2)])
          (cond
            [(fx<? cls1 cls2)
              (cons cls1 (union1d (cdr ls1) ls2))]
            [(fx=? cls1 cls2)
              (cons cls1 (union1d (cdr ls1) (cdr ls2)))]
            [else (cons cls2 (union1d ls1 (cdr ls2)))] )) ])]
    [ (ls1 ls2 ls3)
      (cond
        [(null? ls1) (union1d ls2 ls3) ]
        [(null? ls2) (union1d ls1 ls3) ]
        [(null? ls3) (union1d ls1 ls2) ]
        [else
          (let ([min (fxmin (car ls1) (car ls2) (car ls3))])
            (cons min (union1d (if (fx=? min (car ls1)) (cdr ls1) ls1)
                               (if (fx=? min (car ls2)) (cdr ls2) ls2)
                               (if (fx=? min (car ls3)) (cdr ls3) ls3)))) ]) ]))
                                                                                            ;
(define (union1d! x1s x2s)               ;; {Fixnum}! {Fixnum}! -> {Fixnum}
  ;; produce union1d reusing pairs from x1s and x2s
  (example: (union1d! (list 0 2 4) (list 1 2 3)) => '(0 1 2 3 4) )
  (cond
    [(null? x1s)  x2s ]
    [(null? x2s)  x1s ]
    [else (let ([x1 (car x1s)] [x2 (car x2s)])
            (letrec (
                [out  (if (fx<=? x1 x2)  x1s  x2s ) ]
                [run1 (lambda (x1s x2s tail)       ;; tail is prev x1
                        (cond
                          [(null? x1s)  (set-cdr! tail x2s)  out ]
                          [(null? x2s)  out ]
                          [else (let ([x1 (car x1s)] [x2 (car x2s)])
                              (cond
                                [(fx<? x1 x2)
                                  (run1 (cdr x1s) x2s x1s) ]
                                [(fx>? x1 x2)
                                  (set-cdr! tail x2s)
                                  (run2 x1s (cdr x2s) x2s) ]
                                [else (run1 (cdr x1s) (cdr x2s) x1s) ])) ])) ]
                [run2 (lambda (x1s x2s tail)       ;; tail is prev x2
                        (cond
                          [(null? x2s)  (set-cdr! tail x1s)  out ]
                          [else (let ([x1 (car x1s)] [x2 (car x2s)])
                              (cond
                                [(fx<? x1 x2)
                                  (set-cdr! tail x1s)
                                  (run1 (cdr x1s) x2s x1s) ]
                                [(fx>? x1 x2)
                                  (run2 x1s (cdr x2s) x2s) ]
                                [else (run1 (cdr x1s) (cdr x2s) x1s) ])) ])) ])
              (cond
                [(fx<? x1 x2)  (run1 (cdr x1s) x2s x1s) ]
                [(fx>? x1 x2)  (run2 x1s (cdr x2s) x2s) ]
                [else  (run1 (cdr x1s) (cdr x2s) x1s) ]) )) ]))
                                                                                            ;
(define (setdiff1d l1 l2)                ;; {Fixnum} {Fixnum} -> {Fixnum}
  ;; produce difference of sorted lists of fixnums
  (example: (setdiff1d '(1 2 3 4) '(1 3 5)) => '(2 4) )
  (let loop ([l1 l1] [l2 l2] [result (list)])
    (cond [(null? l1)  (!reverse! result)]
          [(or (null? l2) (fx<? (car l1) (car l2)))
              (loop (cdr l1) l2 (cons (car l1) result))]
          [(fx>? (car l1) (car l2))
              (loop l1 (cdr l2) result)]
          [else (loop (cdr l1) (cdr l2) result)])))
                                                                                            ;
(define (setdiff1d! l1 l2)               ;; {Fixnum}! {Fixnum} -> {Fixnum}
  ;; produce l1 with elements in l2 removed, mutating l1
  (example: (setdiff1d! (list 1 2 3 4) '(1 3 5)) => '(2 4) )
  (example: (setdiff1d! (list 1 2 3 4) '(2 3)  ) => '(1 4) )
  (let skip ([l1 l1] [l2 l2])
    (cond
      [(or (null? l1) (null? l2))  l1 ]
      [(fx=? (car l1) (car l2))  (skip (cdr l1) (cdr l2)) ]
      [else
        (let loop ([c1 l1] [c2 l2] [prev l1])
          (cond
            [(or (null? c1) (null? c2))  l1 ]
            [(fx<? (car c1) (car c2))  (loop (cdr c1) c2 c1) ]
            [(fx>? (car c1) (car c2))  (loop c1 (cdr c2) prev) ]
            [else
              (set-cdr! prev (cdr c1))
              (loop (cdr c1) (cdr c2) prev) ])) ])))
                                                                                            ;
(define (in1d x1s x2s)                   ;; {X} {X} -> Bits
  ;; produce index mask of x1s that are also in x2s
  (example: (in1d '(1 2 8 9 5) '(9 8 1)) => #b01101 )
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
  (example: (include-by-mask '(1 2 3 4) #b1101) => '(1 3 4) )
  (vector->list
    (vector-refs 
      (list->vector xs) 
      (list->vector (bitwise->list mask)))))
                                                                                            ;
(define (exclude-by-mask xs mask)        ;; {X} Bits -> {X}
  ;; omit elements of xs corresponding to 1 bits in mask
  (example: (exclude-by-mask '(1 2 3 4) #b1101) => '(2) )
  (cond
    [(zero? mask) xs]
    [(odd?  mask)
      (exclude-by-mask (cdr xs) (bitwise-arithmetic-shift-right mask 1))]
    [else
      (cons (car xs) (exclude-by-mask (cdr xs) (bitwise-arithmetic-shift-right mask 1)))]))
                                                                                            ;
  (define cores 10)
                                                                                            ;                                                                                            ;
  (define thread-if 3)                   ;; fork if vector-length >= thread-if
                                                                                            ;                                                                                            ;
  (define random-seed-limit  (- (expt 2 32) 1))
                                                                                            ;                                                                                            ;
  (define random-state  (make-thread-parameter 0))
                                                                                            ;                                                                                            ;
(define (new-thread-vector-for-each f . vs);; (X ... -> ) (Vectorof X) ... ->
  ;; in a new thread for each, apply f to elements of vs, return when all finished
  (let* ( [todo         (vector-length (car vs))]
          [max-threads  (fx* 2 cores)]
          [thread-limit (if (fx<=? todo max-threads)  max-threads
                            (fx1+ (fxdiv todo (fx1+ (fxdiv todo max-threads)))))])
    (if (fx<? todo thread-if)
      (apply vector-for-each (lambda xs  ;; could be (apply vector-for-each f vs)
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
            (random-state (fx1+ (random random-seed-limit)))
            (fork-thread (lambda ()                  ;; child thread
                    (random-seed (random-state))     ;; new random-seed for each thread
                    (apply f xs)                     ;;
                    (with-mutex mutex                ;;
                      (set! todo (fx1- todo))        ;;
                      (when (fxzero? todo)           ;;
                        (condition-signal done))     ;;
                      (set! threads (fx1- threads))  ;;
                      (condition-signal free)))      ;;
              ))
          vs)
        (with-mutex mutex
          (unless (fxzero? todo)
            (condition-wait done mutex)))))))
                                                                                            ;
(define-record-type tp                   ;; ThreadPool
  (fields
    (immutable queue)                    ;; (Vector [head|tail]->WorkItem)
    (mutable head)                       ;; [WorkItem is (Proc . Arg),
    (mutable tail)                       ;;  Pair arg is treated as arg-list]
    (mutable pending)                    ;; Nat [outstanding work]
    (immutable lock)                     ;; Mutex
    (immutable work)                     ;; Condition [item in queue]
    (immutable room)                     ;; Condition [room in queue]
    (immutable done))                    ;; Condition [no outstanding items]
  (protocol
    (lambda (new)
      (lambda (bound)
        (new (make-vector bound) 0 0 0 (make-mutex)
          (make-condition) (make-condition) (make-condition))))))
                                                                                            ;
(define (tp-incr tp index)               ;; ThreadPool Nat -> Nat
  ;; produce incremented index in tp's queue
  (fxmod (fx1+ index) (vector-length (tp-queue tp))))
                                                                                            ;
(define (tp-enqueue tp item)             ;; ThreadPool WorkItem ->
  ;; add item to tp's queue (may wait til room)
  (mutex-acquire (tp-lock tp))
  (let loop ()
    (let* ( [tail  (tp-tail tp)]
            [tail^ (tp-incr tp tail)])
      (cond
        [(fx=? tail^ (tp-head tp))
         (condition-wait (tp-room tp) (tp-lock tp))
         (loop) ]
        [else
         (vector-set! (tp-queue tp) tail item)
         (tp-pending-set! tp (fx1+ (tp-pending tp)))
         (tp-tail-set! tp tail^)
         (condition-signal (tp-work tp))
         (mutex-release (tp-lock tp)) ]))))
                                                                                            ;
(define (tp-dequeue tp)                  ;; ThreadPool -> WorkItem
  ;; produce next item from tp's queue (may wait til available)
  (mutex-acquire (tp-lock tp))
  (let loop ()
    (let ([head (tp-head tp)])
      (cond
        [(fx=? head (tp-tail tp))
         (condition-wait (tp-work tp) (tp-lock tp))
         (loop)]
        [else
         (tp-head-set! tp (tp-incr tp head))
         (let ([item (vector-ref (tp-queue tp) head)])
           (condition-signal (tp-room tp))
           (mutex-release (tp-lock tp))
           item) ]))))
                                                                                            ;
(define (tp-wait-til-done tp)            ;; ThreadPool ->
  ;; return when nothing pending in tp's queue
  (mutex-acquire (tp-lock tp))
  (unless (fxzero? (tp-pending tp))
    (condition-wait (tp-done tp) (tp-lock tp)))
  (mutex-release (tp-lock tp)))
                                                                                            ;
(define (make-thread-pool nthreads)      ;; Nat -> ThreadPool
  ;; produce tp with nthreads threads to run (f . arg) from queue
  (let ([tp (make-tp nthreads)])
    (do ([n nthreads (fx1- n)])
        ((fxzero? n))
      (random-state (fx1+ (random random-seed-limit)))
      (fork-thread (lambda ()
          (random-seed (random-state))   ;; new random-seed for each thread
          (let loop ()                   ;; work threads loop indefinitely
            (let* ( [work (tp-dequeue tp)]
                    [arg  (cdr work)])
                (if (pair? arg)
                  (apply (car work) arg)
                  ((car work) arg)))
            (mutex-acquire (tp-lock tp))
            (tp-pending-set! tp (fx1- (tp-pending tp)))
            (when (fxzero? (tp-pending tp))
              (condition-signal (tp-done tp)))
            (mutex-release (tp-lock tp))
            (loop)))))
    tp))
                                                                                            ;
(define (threaded-for-each-x tp f n)     ;; ThreadPool (Nat -> ) Nat ->
  ;; queue each (f i) for i [0..n-1] to tp, return when all done
  (if (fx=? n 1)  (f 0)
      (do ([i 0 (fx1+ i)])
          ((fx=? i n) (tp-wait-til-done tp))
        (tp-enqueue tp (cons f i)))))
                                                                                            ;
#;(define (threaded-for-each-x tp f n)   ;; ThreadPool (Nat -> ) Nat ->
  ;; (unthreaded)
  (do ([i 0 (fx1+ i)])
          ((fx=? i n))
        (f i)))
                                                                                            ;
(define (vectors-ref vs n)               ;; {VectorOfX} Nat -> {X}
  ;; produce nth element of each v
  (map (lambda (v)
      (vector-ref v n))
    vs))
                                                                                            ;
(define (threaded-vector-for-each tp f . vs) ;; ThreadPool ({X} -> ) {(Vectorof X)} ->
  ;; queue each (f . corresponding elements of vs) to tp, return when all done
  (let ([n (vector-length (car vs))])
    (threaded-for-each-x tp
      (lambda (x) (apply f (vectors-ref vs x))) n)))
                                                                                          ;
(define cc-invocations 0)                ;; Fixnum 
(define cc-counter-n 0)
(define cost-center-0 (make-cost-center))
(define cost-center-1 (make-cost-center))
                                                                                            ;
(define (overhead-thunk)                 ;; ->
  ;; count in a thunk to get timing overhead (not VR because they are thread local)
  (set! cc-invocations (fx1+ cc-invocations)))
                                                                                            ;
(define with-thunk-time                  ;; ( -> ? ) [ ( -> Fixnum) ] -> ?
  ;; track thunk cpu, cc overhead (optimization may make overhead offset inaccurate)
  (case-lambda
    [(thunk)
      (with-thunk-time thunk (lambda () 1)) ]
    [(thunk count)
      (with-cost-center #t cost-center-0 overhead-thunk)
      (let ([result (with-cost-center #t cost-center-1 thunk)])
        (set! cc-counter-n (fx+ cc-counter-n (count)))
        result)  ]))
                                                                                            ;
(define (reset-thunk-time)               ;; ->
  (set! cc-invocations 0)
  (set! cc-counter-n 0)
  (reset-cost-center! cost-center-0)
  (reset-cost-center! cost-center-1))
                                                                                            ;
(define (display-thunk-time)             ;; ->
  (when (fxpositive? cc-invocations)
    (let* ( [t0    (cost-center-time cost-center-0)]
            [t1    (cost-center-time cost-center-1)]
            [t     (time-difference t1 t0)]
            [secs  (time-second t)]
            [ns    (time-nanosecond t)]
            [tns   (fx+ (fx* secs 1000000000) ns)]
            [t/x   (fxdiv tns cc-invocations)]
            [alloc (cost-center-allocation-count cost-center-1)]
            [allst (number->string alloc)]
            [nch   (string-length allst)]
            [bytes (substring allst (max 0 (- nch 6)) nch)]
            [mb    (substring allst 0 (max 0 (- nch 6)))])
        (for-each display `(
            ,secs "." ,(fxdiv ns 10000000) "s "
            "for " ,cc-invocations " xs ("
              ,(if (fx<? t/x 1000)
                (string-append (number->string t/x) "ns/x")
                (string-append
                  (number->string (/ (fxdiv (fx+ t/x 50) 100) 10.)) "µs/x"))
              ,(if (zero? cc-counter-n)  ""
                (string-append "; " (number->string cc-counter-n)))
              ,(if (zero? alloc)  ""
                (string-append "; " mb " " bytes))
              "); ")))
    (when (fxpositive? (vr0@))
      (for-each display `(
          ,(vr0@) " " ,(vr1@) "; ")))
    (newline)))
                                                                                            ;
(define (vr0@)   (virtual-register 0))
(define (vr1@)   (virtual-register 1))
(define (vr0! x) (set-virtual-register! 0 x))
(define (vr1! x) (set-virtual-register! 1 x))
                                                                                            ;
  (define first-n 0)
                                                                                            ;
  (define first-mutex (make-mutex))
                                                                                            ;
(define (first n thunk)                  ;; Fixnum Thunk ->
  ;; call thunk first n times (for debug output)
  (when (fx<? first-n n)
    (set! first-n (fx1+ first-n))
    (with-mutex first-mutex
      (thunk))))

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
                                                                                            ;
(define expect-count (vector 0))
                                                                                            ;
(define (expects)
  ;; produce count of expects, reset count
  (let ([n (vector-ref expect-count 0)])
    (vector-set! expect-count 0 0)
    n))
                                                                                            ;
(define-syntax expect                    ;; ((X ... -> Y) X ...) Y -> Error |
  ;; check that function application(s) to arguments match expected values
  (lambda (x)                            
    (syntax-case x ()                    ;; [expect ([fn args] expected ) ... ]
      [ (_ [expr expected] ...)          ;; expr matches [fn args]
        #'(begin
            (let ([result expr])           ;; eval expr once, no output if check passes
              #;
              (when (equal? result expected) (display "."))
              (unless (equal? result expected)
                (for-each display
                  `("\n**" expr #\newline
                    "  expected: " ,expected #\newline
                    "  returned: " ,result  #\newline))
                (exit))
              (vector-set! expect-count 0 (fx1+ (vector-ref expect-count 0)))) ...) ]
      )))

;; --- Smoke tests ---
                                                                                            ;
  (let* ( [source (build-vector 100 id)]
          [select (vector->list [vector-sample source 50])])
    [expect ( [length select] 50 )
            ( [list? (for-all memq select (make-list 50 (vector->list source)))] #t )  ;; all from source
            ( [map (lambda (x) (length (remq x select))) select] (make-list 50 49)) ]) ;; all different
  [expect ( [(lambda (x . args)          ;; example fn with kws k1 & k2: default is (x 1 2)
              (append (list x) [key-word-args args '([k1 . 1] [k2 . 2])] ))
            99 '[k2 . 22] '[k1 . 11] '[k3 . 33] ] '(99 11 22) )]           ;; apply fn to 99 '[k2 ...  
  #;
  (let ([xs (make-vector 3 #f)])
    [expect ( [begin [threaded-vector-for-each (lambda (x ix)
                          (vector-set! xs ix (+ x 1)))
                        '#(1 2 3) '#(0 1 2)]
                      xs]   '#(2 3 4))])
  
(vr0! 0)
(vr1! 0)

;; (not effective on import, but will run when any export used):
;; (eval-when (compile eval load visit revisit)
     (check-examples)
;; )

)

#| Notices:

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <https://www.gnu.org/licenses/>.
  
  License: <https://www.gnu.org/licenses/agpl-3.0.txt>
  Contact: <https://github.com/rogerturner/HTM-scheme/issues/new/choose>  |#
