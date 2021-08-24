;; === HTM-scheme Coordinates  (C) 2021 Roger Turner. ===
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
                        _______
            ^ r        / . . . \         3 hexagonal "cortical columns"
           /          / . . . . \        (19 minicolumns each - would normally be 127)
          /    ______/ . . . . . )
         /    / . . . \ . . . . /        triangular minicolumn lattice is actually
        /    / . . . . \ . . . /         perfectly regular (no extra gap between ccs)
       /    ( . . . . . ) . . . \
      /      \ . . . . / . . . . \
     /        \ . . . / . . . . . )
    /          ------- \ . . . . /
   /                    \ . . . /                   
  /                      -------
 o---------------------------------> q

"Axial" axes: q is like x axis, r replaces y axis but is at 60º to q
In axial coords, distance squared from origin to (q,r) is q^2 + r^2 + q*r

See htm_concept.ss for type and data structure description and code conventions.

  |#
  
  #!chezscheme

(library (HTM-scheme HTM-scheme math coordinates)
                                                                                            ;
(export
q-coord-of-cc-centre
r-coord-of-cc-centre
q-coord-of-minicol
r-coord-of-minicol
within-cc-distance2
cc-distance2
mc-distance2
  )
                                                                                            ;
(import
  (except (chezscheme) reset)
  (parameters)
  (HTM-scheme HTM-scheme algorithms htm_prelude)
  (HTM-scheme HTM-scheme algorithms htm_concept))
                                                                                            ;
  (implicit-exports #f)

;; === Hexagonal coordinate system ===
                                                                                            ;
(define (hex-lattice radius tile)        ;; Nat Nat -> {(Fixnum . Fixnum)}
  ;; Produce hexagonal lattice of (q,r) coordinate points with given radius
  ;; Head of result list is origin, points follow in sequence of hexagons
  ;; of increasing size (last 6*radius points are boundary of lattice)
  ;; radii 0, 1, 2, 3, ... produce 1, 7, 19, 37 etc points
  ;; positive tile => lattice of centres of sub-hexagons with that radius
  (define (hexagon)                      ;; -> {(q . r)}
    (define (side q+ r+ base)            ;; Nat Nat (q . r) -> {(q . r)}
      (define (move-from prev)           ;; (q . r) -> (q . r)
        ;; produce next point in direction (q+,r+) from previous point
        (cons (+ q+ (car prev)) (+ r+ (cdr prev))))
      ;; (side) produce side in direction (q+,r+) from base
      (do ( [i 1 (+ i 1)]
            [side (list (move-from base)) (cons (move-from (car side)) side)])
          ((= i radius) side)))
    (define (add-side q+ r+ prev-side)   ;; Nat Nat {(q . r)} -> {(q . r)}
      ;; produce next side using last point of previous side as base
      (append (side q+ r+ (car prev-side)) prev-side))
    ;; (hexagon) produce points of one hexagon using radius and tile
    (let* (
        [start  (cons (* radius (+ tile 1)) (* radius tile))]
        [points (side     (- (+ (* 2 tile) 1)) (+ tile 1) start)]
        [points (add-side (- (+ tile 1))       (- tile)              points)]
        [points (add-side tile                 (- (+ (* 2 tile) 1))  points)]
        [points (add-side (+ (* 2 tile) 1)     (- (+ tile 1))        points)]
        [points (add-side (+ tile 1)           tile                  points)]
        [points (add-side (- tile)             (+ (* 2 tile) 1)      points)])
      points))
  ;; (hex-lattice) produce lattice by nesting smaller lattice in outer hexagon
  (if (zero? radius)  (list '(0 . 0))
      (append (hex-lattice (- radius 1) tile) (hexagon))))
                                                                                            ;
  ;; local q,r coords for mc within cc tile; number of minicolumns <= 397 (radius 11)
  (define mc-q (sint-list->bytevector (map car (hex-lattice 11 0)) (native-endianness) 1))
  (define mc-r (sint-list->bytevector (map cdr (hex-lattice 11 0)) (native-endianness) 1))
                                                                                            ;
  (alias mcref bytevector-s8-ref)
                                                                                            ;
  (define cc-radius                      ;; Nat: 
    (case minicolumns/macrocolumn
      [( 91)  5 ]
      [(127)  6 ]
      [(169)  7 ]
      [(217)  8 ]
      [(271)  9 ]
      [(331) 10 ]
      [(397) 11 ]
      [(469) 12 ]
      [(547) 13 ]
      [(631) 14 ]
      [(721) 15 ]
      [(817) 16 ]
      [(919) 17 ]
      [(1027) 18 ]
      [else   0 ]))
                                                                                            ;
  (define bytes
    (if (< minicolumns/macrocolumn 397) 1 2))
                                                                                            ;
  ;; cc centre coords for lattice of up to 127 ccs of size minicolumns/cc                                                                                          
  (define cc-q                           ;; Bytevector: cortical column centre q coords
    (sint-list->bytevector
      (map car (hex-lattice 6 cc-radius))
      (native-endianness) bytes))
                                                                                            ;
  (define cc-r                           ;; Bytevector: cortical column centre r coords
    (sint-list->bytevector
      (map cdr (hex-lattice 6 cc-radius))
      (native-endianness) bytes))
                                                                                            ;
  (define-syntax ccref
    (lambda (x)
      (syntax-case x ()
        [ (_ bv bvx)
          (if (< minicolumns/macrocolumn 397)
            #'(bytevector-s8-ref bv bvx)
            #'(bytevector-s16-native-ref bv (fx* 2 bvx))) ])))
                                                                                            ;
(define (q-coord-of-cc-centre ccx)       ;; CCX -> Fixnum
  (ccref cc-q ccx))
                                                                                            ;
(define (r-coord-of-cc-centre ccx)       ;; CCX -> Fixnum
  (ccref cc-r ccx))
                                                                                            ;
(define (q-coord-of-minicol mcx)         ;; ColX -> Fixnum
  (mcref mc-q mcx))
                                                                                            ;
(define (r-coord-of-minicol mcx)         ;; ColX -> Fixnum
  (mcref mc-r mcx))
                                                                                            ;
(define (distance2 q-diff r-diff)        ;; Fixnum Fixnum -> Nat
  ;; produce distance squared between points with coord differences (q-diff, r-diff)
  (fx+ (fx* q-diff (fx+ q-diff r-diff)) (fx* r-diff r-diff)))
                                                                                            ;
(define (within-cc-distance2 mx1 mx2)    ;; ColX ColX -> Nat
  ;; produce distance squared between minicolumns mx1 and mx2
  (let ([q-diff (fx- (mcref mc-q mx1) (mcref mc-q mx2))]
        [r-diff (fx- (mcref mc-r mx1) (mcref mc-r mx2))])
    (distance2 q-diff r-diff)))
                                                                                            ;
(define (cc-distance2 cx1 cx2)           ;; CCX CCX -> Nat
  ;; produce distance squared between centres of cortical columns cx1 and cx2
  (let ([q-diff (fx- (ccref cc-q cx1) (ccref cc-q cx2))]
        [r-diff (fx- (ccref cc-r cx1) (ccref cc-r cx2))])
    (distance2 q-diff r-diff)))
                                                                                            ;
(define (mc-distance2 cx1 mx1 cx2 mx2)   ;; CCX ColX CCX ColX -> Nat
  ;; produce distance squared between minicolumns (mx) within cortical column (cx)
  (let ([mc-diff-q (fx- (mcref mc-q mx1) (mcref mc-q mx2))]
        [mc-diff-r (fx- (mcref mc-r mx1) (mcref mc-r mx2))]
        [cc-diff-q (fx- (ccref cc-q cx1) (ccref cc-q cx2))]
        [cc-diff-r (fx- (ccref cc-r cx1) (ccref cc-r cx2))])
    (let ([q-diff (fx+ mc-diff-q cc-diff-q)]
          [r-diff (fx+ mc-diff-r cc-diff-r)])
      (distance2 q-diff r-diff))))

;; -- Smoke tests --
                                                                                            ;
  [expect [(hex-lattice 0 0)  '((0 . 0)) ]]
  [expect [(hex-lattice 1 0)  '((0 . 0)
                                (1 . 0) (1 . -1) (0 . -1) (-1 . 0) (-1 . 1) (0 . 1)) ]]
  [expect [(take 8 (hex-lattice 2 0))
                              '((0 . 0)
                                (1 . 0) (1 . -1) (0 . -1) (-1 . 0) (-1 . 1) (0 . 1)
                                (2 . 0)) ]]
  [expect [(list-tail (hex-lattice 2 0) 17)
                              '((0 . 2) (1 . 1)) ]]
  [expect [(hex-lattice 1 1)  '((0 . 0)
                                (2 . 1) (3 . -2) (1 . -3) (-2 . -1) (-3 . 2) (-1 . 3)) ]]
  [expect [(hex-lattice 1 2)  '((0 . 0)
                                (3 . 2) (5 . -3) (2 . -5) (-3 . -2) (-5 . 3) (-2 . 5)) ]]

)