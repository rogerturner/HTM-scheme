#!chezscheme

;; ====== HTM-scheme Coordinate Encoder Copyright 2020 Roger Turner. ======
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
  #|

Given an integer coordinate in a 2-dimensional space, and a radius around
that coordinate, returns an SDR representation of that position.
Translated from NuPIC coordinate.py, see comments there for more info.
  Description of algorithm from .py:
  1. Find all the coordinates around the input coordinate, within the
     specified radius.
  2. For each coordinate, use a uniform hash function to
     deterministically map it to an "order".
  3. Of these coordinates, pick the top W by order, where W is the
     number of active bits desired in the SDR.
  4. For each of these W coordinates, use a uniform hash function to
     deterministically map it to one of the bits in the SDR. Make this bit
     active.
  5. This results in a final SDR with exactly W bits active.
  |#

(library (HTM-scheme HTM-scheme encoders coordinate)
                                                                                            ;
(export
  make-ce
  ce-w
  ce-w-set!
  ce-n
  ce-n-set!
  ce-radius
  ce-radius-set!
  encode)
                                                                                            ;
(import 
  (chezscheme)
  (HTM-scheme HTM-scheme algorithms htm_prelude)
  (HTM-scheme HTM-scheme algorithms htm_concept))

;; Types
;; Coord    Scheme Complex with Fixnum real and imag parts

(define-record-type ce                   ;; CoordinateEncoder [CE]
  (fields
    (mutable w)                          ;; Nat: number of 1 bits in SDR
    (mutable n)                          ;; Nat: size of SDR
    (mutable radius)                     ;; Nat: 
    )
  (sealed #t) (opaque #t) (nongenerative ce))

(define (encode ce centre)               ;; CE Coord -> SDR
  ;; produce SDR for centre using ce
#;> (define (_order-for-coord c)         ;; Coord -> Nat
      ;; produce order for c 
      (_hash-coordinate c)
      #;(let* ([xd (fx- (real-part c) (real-part centre))] ;; gaussian weight by distance from centre?
             [yd (fx- (imag-part c) (imag-part centre))]
             [d  (exact (round (sqrt (fx+ (fx* xd xd) (fx* yd yd)))))]
             [m  (if (fx>? d 9) 1
                   (case d
                     [(9) 2]
                     [(8) 4]
                     [(7) 7]
                     [(6) 10]
                     [(5) 14]
                     [(4) 18]
                     [(3) 21]
                     [(2) 23]
                     [(1) 24]
                     [(0) 25]))])
        (fx* m (_hash-coordinate c))))
#;> (define (_top-coords neighborhood w) ;; (Vectorof Coord) Nat -> (Vectorof Coord)
      ;; produce w top elements of neighborhood by order
      (let* ( [orders  (vector-map _order-for-coord neighborhood)]
              [indices (vector-take w (vector-sort (lambda (o1 o2)
                            (fx>? (vector-ref orders o1) (vector-ref orders o2)))
                          (indexes neighborhood)))])
      (vector-refs neighborhood indices)))
  (let* ( [neighbors (_neighbors centre (ce-radius ce))]
          [winners   (_top-coords neighbors (ce-w ce))]
          [indices   (vector-map (lambda (c)
                          (_bit-for-coord c (ce-n ce)))
                        winners)])
    (vector->bitwise indices)))
                                                                                            ;
(define (_neighbors centre radius)       ;; Coord Nat -> (Vectorof Coord)
  ;; produce coordinates within radius of centre (rounded for radius >= 5)
  (let* (
      [width  (fx1+ (fx* 2 radius))]
      [result (make-vector (fx- (fx* width width)
                                (cond [(fx<? radius 5) 0 ]
                                      [(fx<? radius 7) 32]
                                      [else 68])))])
    (let ((xc (real-part centre))
          (yc (imag-part centre)))
      (let xloop ([x (- radius)] [ri 0])
        (let* (
            [corner (if (fx<? radius 5)  4
                        (fx- radius (fxabs x)))]
            [corner (if (fx<? radius 7)
                (case corner [(0) 4] [(1) 2] [(2 3) 1] [else 0])
                (case corner [(0) 6] [(1) 4] [(2) 3] [(3) 2] [(4 5) 1] [else 0]))])
          (if (fx<=? x radius)
            (xloop (fx1+ x)
              (let yloop ([y (fx- corner radius)]
                          [ri ri])
                (cond
                  [ (fx<=? y (fx- radius corner))
                      (vector-set! result ri 
                        (make-rectangular (fx+ xc x) (fx+ yc y)))
                      (yloop (fx1+ y) (fx1+ ri))]
                  [ else ri])))
            result))))))
                                                                                            ;
(define hashes                           ;; (Vectorof Nat)
  ;; random permutation of 0..256k
  (let* ( [k256 (* 512 512)]
          [v    (vector-sample (build-vector k256 id) k256)]
          [h    (make-bytevector (* 4 k256))])
    (do ((i (fx1- k256) (fx1- i))) ((fxnegative? i) h)
      (bytevector-u32-native-set! h (fx* i 4) (vector-ref v i)))))
                                                                                            ;
(define (spread9 x)                      ;; Fixnum -> Fixnum
  ;; produce 18 bits from 9 low-order bits of x: bit 0->0, 1->2, ...8->16
  (let* (
      [x (fxior (fxarithmetic-shift-left (fxand x #xFF00) 8) (fxand x #x00FF))]
      [x (fxand (fxior (fxarithmetic-shift-left x 4) x) #x0F0F0F0F)]
      [x (fxand (fxior (fxarithmetic-shift-left x 2) x) #x33333333)]
      [x (fxand (fxior (fxarithmetic-shift-left x 1) x) #x55555555)])
    x ))
                                                                                            ;
(define (_hash-coordinate c)             ;; Coord -> Nat
  ;; produce uniformly distributed (?) hash of c (9 bit coords)
  (bytevector-u32-native-ref hashes 
      (fx* 4 (fxand #x3FFFF (fxior
          (spread9 (real-part c))
          (fxarithmetic-shift-left (spread9 (imag-part c)) 1))))))
                                                                                            ;
(define (_bit-for-coord c n)             ;; Coord Nat -> Nat
  ;; produce value in [0..n> for c
  (random-seed (fx* 997 (_hash-coordinate c)))
  (random n))

;; Smoke tests
                                                                                            ;
  (expect ([_neighbors (make-rectangular 2 2) 1]
            '#(1+1i 1+2i 1+3i 2+1i 2+2i 2+3i 3+1i 3+2i 3+3i)))
  
)