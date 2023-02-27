;; © 2022 Roger Turner <https://github.com/rogerturner/HTM-scheme/issues/new/choose>
;; SPDX-License-Identifier: AGPL-3.0-or-later  (see Notices below)

#| HTM-scheme Sorting (Customized sort functions)

Uniquing in-place merge sorts for Fixnums (SDRs) and Segments.
Derived from [Shivers 2002 A simple and efficient natural merge sort, Appendix D]
and small-list sorting from chezscheme library.

  |#
  
  #!chezscheme

(library (frameworks sorting)

(export
sort-unique!
sort-unique-by!
  )

(import
  (chezscheme)
  (frameworks htm-prelude))
                                                                                            ;
  (implicit-exports #f)

(define (fxgetrun1! fxs)                 ;; {Fixnum}! -> Fixnum Pair {Fixnum}
  ;; produce run length, last pair, following list; omit duplicates
  (example:
    (let*-values (
      [(fxs)        (list 4 5 5 1 2)]
      [(*5 *1)      (values (cdr fxs) (cdddr fxs))]
      [(rl lp rest) (fxgetrun1! fxs)])
        (and (= rl 2) (eq? lp *5) (eq? rest *1) (equal? fxs '(4 5 5 1 2)))) )
  (example:
    (let*-values (
      [(fxs)        (list 4 5 5 6 1 2)]
      [(*6 *1)      (values (cdddr fxs) (cddddr fxs))]
      [(rl lp rest) (fxgetrun1! fxs)])
        (and (= rl 3) (eq? lp *6) (eq? rest *1) (equal? fxs '(4 5 6 1 2)))) )
  (let ->run ([fxs fxs] [x (car fxs)] [i 1] [next (cdr fxs)])
    (if (pair? next)
      (let ([y (car next)])
        (cond
          [(fx<? y x) (values i fxs next) ]
          [(fx=? y x) ;; (set-cdr! fxs ...)
            (let ->skip ([next (cdr next)])
              (if (pair? next)
                (let ([y (car next)])
                  (cond
                    [(fx=? y x) (->skip (cdr next)) ]
                    [(fx<? y x) (values i fxs next) ]
                    [else
                      (set-cdr! fxs next)
                      (->run next y (fx1+ i) (cdr next)) ]))
                (values i fxs next))) ]
          [else (->run next y (fx1+ i) (cdr next)) ]))
      (values i fxs next))))
                                                                                            ;
(define (fxmerge1! a enda b endb)        ;; {Fixnum} Pair {Fixnum} Pair -> {Fixnum} Pair
  ;; produce merge of a, b where enda, endb are their last pairs; omit duplicates
  (example:
    (let*-values (
        [(a    b)    (values (list 1 3) (list 2 4))]
        [(enda endb) (values (cdr a) (cdr b))]
        [(ab endab)  (fxmerge1! a enda b endb)])
      (and (eq? ab a) (eq? endab endb) (equal? ab '(1 2 3 4)))))
  (example:
    (let*-values (
        [(a    b)    (values (list 1 3) (list 1 4))]
        [(enda endb) (values (cdr a) (cdr b))]
        [(ab endab)  (fxmerge1! a enda b endb)])
      (and (eq? ab a) (eq? endab endb) (equal? ab '(1 3 4)))))
  (example:
    (let*-values (
        [(a    b)    (values (list 1 2) (list 2 4))]
        [(enda endb) (values (cdr a) (cdr b))]
        [(ab endab)  (fxmerge1! a enda b endb)])
      (and (eq? ab a) (eq? endab endb) (equal? ab '(1 2 4)))))
  (example:
    (let*-values (
        [(a    b)    (values (list 1 4) (list 2 4))]
        [(enda endb) (values (cdr a) (cdr b))]
        [(ab endab)  (fxmerge1! a enda b endb)])
      (and (eq? ab a) (eq? endab enda) (equal? ab '(1 2 4)))))
  (letrec (
      [scan-a (lambda (prev x a y b)     ;; x is (car a), y is (car b)
          (cond                          ;; zip down a until we..
            [(fx<? y x)                  ;; find an elt > (car b)
              (maybe-set-cdr! prev b)
              (let ([next-b (cdr b)])
                (if (eq? b endb)
                  (begin (maybe-set-cdr! b a) enda)      ;; done
                  (scan-b b x a (car next-b) next-b))) ]
            [(fx=? y x)
              (if (eq? a enda)
                (begin (maybe-set-cdr! prev b) endb)
                (let ([next-a (cdr a)])
                  (scan-a prev (car next-a) next-a y b))) ]
            [(eq? a enda)
              (maybe-set-cdr! a b) endb ]                ;; done
            [else
              (let ([next-a (cdr a)])    ;; continue scan
                (scan-a a (car next-a) next-a y b)) ]))]
                                                                                            ;
      [scan-b (lambda (prev x a y b)
          (cond                          ;; zip down b while its..
            [(fx<? y x)                  ;; elts are < (car a)
              (if (eq? b endb)
                (begin (maybe-set-cdr! b a) enda)        ;; done
                (let ([next-b (cdr b)])  ;; continue scan
                  (scan-b b x a (car next-b) next-b))) ]
            [(fx=? y x)
              (if (eq? b endb)
                (begin (maybe-set-cdr! prev a) enda)
                (let ([next-b (cdr b)])
                  (scan-b prev x a (car next-b) next-b))) ]
            [else
              (set-cdr! prev a)
              (if (eq? a enda)
                (begin (maybe-set-cdr! a b) endb)        ;; done
                (let ([next-a (cdr a)])
                  (scan-a a (car next-a) next-a y b))) ])) ]
                                                                                            ;
      [maybe-set-cdr! (lambda (pair val)
          (unless (eq? (cdr pair) val)
            (set-cdr! pair val))) ])
    ;; (fxmerge1! a enda b endb)
    (let ([x (car a)] [y (car b)])
      (cond
        [(fx<? y x)
          (values b                      ;; b starts the result
            (if (eq? b endb)
              (begin (maybe-set-cdr! b a) enda)
              (let ([next-b (cdr b)])
                (scan-b b x a (car next-b) next-b)))) ]
        [(fx=? y x)
          (values a                      ;; a starts the result
            (if (eq? b endb)
              enda
              (let ([next-b (cdr b)])
                (scan-b b x a (car next-b) next-b)))) ]
        [else
          (values a                      ;; a starts the result
            (if (eq? a enda)
              (begin (maybe-set-cdr! a b) endb)
              (let ([next-a (cdr a)])
                (scan-a a (car next-a) next-a y b)))) ]))))
                                                                                            ;
(define (fxgrow1! s ends ls ls2 u lw)    ;; {Fx}! Pair Fx Fx {Fx} Fx -> {Fx} Pair Fx {Fx}
  ;; s is sorted, length ls > 1, with last pair ends
  ;; ls2 is power of 2 <= ls, u is unsorted, lw is positive
  ;; produce sorted result a, enda, la (at least lw), v (unused tail of u)
  (let recur ([s s] [ends ends] [ls ls] [ls2 ls2] [u u] [lw lw])
    (if (and (pair? u) (fx<? ls lw))
      (let*-values (
          [(ls2)  (let ->*2 ([ls2 ls2])
                    (let ([ls2*2 (fx+ ls2 ls2)])
                      (if (fx<=? ls2*2 ls)  (->*2 ls2*2)  ls2)))]
          ;; ls2 is now the largest power of two <= ls
          [(lr endr u2)    (fxgetrun1! u)] ;; u to endr now dup-free
          [(t endt lt u3)  (recur u endr lr 1 u2 ls2)]
          [(st end-st)     (fxmerge1! s ends t endt)])
        (recur st end-st (fx+ ls lt) (fx+ ls2 ls2) u3 lw))
      (values s ends ls u))))
                                                                                            ;
(define (fxsort1! fxs)                   ;; {Fixnum}! -> {Fixnum}
  ;; ascending sort of fxs, without duplicates
  (example:
    (let* ( [fxs (list 1 2 1)]
            [fx^ (fxsort1! fxs)])
      (and (equal? fx^ '(1 2)) (eq? fx^ fxs))))
  (example:
    (let* ( [fxs (list 3 1 2)]
            [*1  (cdr fxs)]
            [fx^ (fxsort1! fxs)])
      (and (equal? fx^ '(1 2 3)) (eq? fx^ *1))))
  (example:
    (let* ( [fxs (list 3 1 1 4 4 2 4 3)]
            [*1  (cdr fxs)]
            [fx^ (fxsort1! fxs)])
      (and (equal? fx^ '(1 2 3 4)) (eq? fx^ *1))))
  (example:
    (let* ( [fxs (list 1 9 0 1 1 1 7 5 3 0 2 4 6 1 8 9 9 9 9)]
            [*0  (cddr fxs)]
            [fx^ (fxsort1! fxs)])
      (and (equal? fx^ '(0 1 2 3 4 5 6 7 8 9)) (eq? fx^ *0))))
  (if (pair? fxs)
    (let*-values (
        [(lr endr rest)  (fxgetrun1! fxs)]  ;; fxs to endr now dup-free
        [(a enda la v)   (fxgrow1! fxs endr lr 1 rest (greatest-fixnum))])
      (set-cdr! enda '())
      a)
    fxs))
                                                                                            ;
(define sort-unique!                     ;; {Fixnum} [Nat] -> {Fixnum}
  ;; ascending sort of fxs [of length n] without duplicates, reusing pairs
  (case-lambda
    [(fxs)
      (fxsort1! fxs) ]
    [(fxs n)
      (fxsort1! fxs) ]))
                                                                                            ;
(define (sort-unique-by! elt< xs)        ;; (X X -> T?) {X} -> {X}
  ;; produce sort of xs by elt<, without duplicates, reusing pairs
  (let ([n (length xs)])
    (if (fx<=? n 1)  xs
        (dosort-by! elt< xs n (cons '() '())))))
                                                                                            ;
(define (dosort-by! elt< xs n loc)       ;; (X X -> T?) {X} Nat Pair -> {X}
  ;; n is length of xs, loc is a list head for domerge-by!
  (cond
    [(fx=? n 2)
      (let* ([xs2 (cdr xs)] [x1 (car xs)] [x2 (car xs2)])
        (cond
          [(eq?  x1 x2) xs2 ]
          [(elt< x1 x2) xs  ]
          [else (set-car! xs x2) (set-car! xs2 x1) xs ])) ]
    [else
      (let* ( [i    (fxsrl n 1)]
              [last (let next ([xs xs] [i i])
                      (if (fx=? i 1)  xs
                          (next (cdr xs) (fx1- i)))) ]
              [rest (cdr last)])
        (set-cdr! last '())
        (domerge-by! elt<
          (if (fx=? i 1)  xs
              (dosort-by! elt< xs i loc))
          (let ([n (fx- n i)])
            (if (fx=? n 1)  rest)
                (dosort-by! elt< rest n loc))
          loc)) ]))
                                                                                            ;
(define (domerge-by! elt< xs1 xs2 loc)   ;; (X X -> T?) {X} {X} Pair -> {X}
  ;; produce merge of xs1 and xs2, omitting duplicates
  (let loop ([xs1 xs1] [xs2 xs2] [loc loc])
    (cond
      [(null? xs1) (set-cdr! loc xs2) ]
      [(null? xs2) (set-cdr! loc xs1) ]
      [else (let ([x1 (car xs1)] [x2 (car xs2)])
        (cond
          [(eq? x2 x1)                   ;; (eq? because used to sort segments)
            (set-cdr! loc xs2)
            (loop (cdr xs1) (cdr xs2) xs2)]
          [(elt< x2 x1)
            (set-cdr! loc xs2)
            (loop xs1       (cdr xs2) xs2)]
          [else
            (set-cdr! loc xs1)
            (loop (cdr xs1) xs2       xs1)])) ]))
  (cdr loc))

;; -- Smoke tests --
                                                                                            ;
  [expect
    ([sort-unique! (list)]         '() )
    ([sort-unique! (list 1 2 3)]    '(1 2 3) )
    ([sort-unique! (list 2 3 1) 3]  '(1 2 3) )
    ([sort-unique! (list 1 2 2 3)]  '(1 2 3) )
    ([sort-unique! (list 0 9 1 1 1 1 7 5 3 0 2 4 6 1 8 9 9 9 9)]  '(0 1 2 3 4 5 6 7 8 9) )]
                                                                                            ;
  (let ()
    (define (s< s1 s2)
      (string<? (symbol->string s1) (symbol->string s2)))
    [expect
      ([sort-unique-by! s< (list)]             '() )
      ([sort-unique-by! s< (list 'a 'b 'c)]    '(a b c) )
      ([sort-unique-by! s< (list 'b 'c 'a)]    '(a b c) )
      ([sort-unique-by! s< (list 'a 'b 'b 'c)] '(a b c) )
      ([sort-unique-by! s< (list 'a 'c 'e 'g 'c 'c 'b 'a 'd 'f 'h 'h)]  '(a b c d e f g h) )] )

;; (not effective on import, but will run when any export used):
;; (eval-when (compile eval load visit revisit)
     (check-examples)
;; )

#| (Earlier sort-unique!, derived from small-list sorting from chezscheme library,
    for Fixnums and omitting duplicates, dolsort! unrolled)
                                                                                            ;
(define sort-unique!                     ;; {Fixnum} [Nat] -> {Fixnum}
  ;; produce ascending sort of ls [of length n] without duplicates, reusing pairs
  (case-lambda
    [(ls)
      (sort-unique! ls (length ls)) ]
    [(ls n)
      (if (fx<=? n 1)  ls
          (dolsort! ls n (cons '() '()))) ]))
                                                                                            ;
  #|  (define (dolsort! elt< ls n loc)   ;; original dolsort! from s/5_6.ss:
        (if (fx= n 1)
            (begin (set-cdr! ls '()) ls)
            (let ([i (fxsrl n 1)])
              (let ([tail (list-tail ls i)])
                (dolmerge! elt<
                  (dolsort! elt< ls i loc)
                  (dolsort! elt< tail (fx- n i) loc)
                  loc)))))  |#
                                                                                            ;
(define (dolsort! ls n loc)              ;; {Fixnum} Nat Pair -> {Fixnum}
  ;; n is length of ls, loc is a list head for dolmerge!
  (cond
    [(fx=? n 2)
      (let* ([ls2 (cdr ls)] [x1 (car ls)] [x2 (car ls2)])
        (cond
          [(fx<? x1 x2) ls  ]
          [(fx=? x1 x2) ls2 ]
          [else (set-car! ls x2) (set-car! ls2 x1) ls ])) ]
    [else
      (let* ( [i    (fxsrl n 1)]
              [last (let next ([ls ls] [i i])
                      (if (fx=? i 1)  ls
                          (next (cdr ls) (fx1- i)))) ]
              [rest (cdr last)])
        (set-cdr! last '())
        (dolmerge!
          (if (fx=? i 1)  ls
              (dolsort! ls i loc))
          (let ([n (fx- n i)])
            (if (fx=? n 1)  rest)
                (dolsort! rest n loc))
          loc)) ]))
                                                                                            ;
(define (dolmerge! ls1 ls2 loc)          ;; {Fixnum} {Fixnum} Pair -> {Fixnum}
  ;; produce merge of ls1 and ls2, omitting duplicates
  (let loop ([ls1 ls1] [ls2 ls2] [loc loc])
    (cond
      [(null? ls1) (set-cdr! loc ls2)]
      [(null? ls2) (set-cdr! loc ls1)]
      [else (let ([cls1 (car ls1)] [cls2 (car ls2)])
        (cond
          [(fx<? cls2 cls1)
            (set-cdr! loc ls2)
            (loop ls1       (cdr ls2) ls2)]
          [(fx=? cls2 cls1)
            (set-cdr! loc ls2)
            (loop (cdr ls1) (cdr ls2) ls2)]
          [else
            (set-cdr! loc ls1)
            (loop (cdr ls1) ls2       ls1)])) ]))
  (cdr loc))
                                                                                            ;
|#

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
