#!chezscheme

;; === HTM-scheme KropffTreves2008 project Copyright 2020 Roger Turner. ===
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Based on code from "ctrl-z-9000-times/KropffTreves2008_reproduction"  ;;
  ;; which is (c) 2018 David McDougall; see license there for permissions. ;;
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
  
Partial replication of experiments reported in:
[McDougall, D 2018 "Reproduction of Kropff & Treves, 2008"]
(https://github.com/ctrl-z-9000-times/KropffTreves2008_reproduction)
  |#

(import
  (except (chezscheme) add1 make-list reset)
  (except (HTM-scheme HTM-scheme algorithms htm_prelude) #;random)
          (HTM-scheme HTM-scheme algorithms htm_concept)
  (prefix (HTM-scheme HTM-scheme encoders   coordinate)                  ce:)
  (prefix (HTM-scheme HTM-scheme algorithms spatial_pooler)              sp:)
  (prefix (HTM-scheme projects KropffTreves2008_reproduction grid_cells) gc:))

;; Types
;; Position   Complex, real and imag parts are Flonums for x, y

  (define pi 3.14159265358979323846264)

(define-record-type env                  ;; Environment [Env]
  ;; square with side |size|,
  ;; tracking |course| of agent at |position| moving with |speed|, |angle|
  (fields
    size                                 ;; Flonum
    (mutable position)                   ;; Position
    (mutable speed)                      ;; Flonum
    (mutable angle)                      ;; Flonum
    (mutable course)                     ;; (Listof Position)
    (mutable seed)
    )
  (protocol
    (lambda (new)
      (lambda (size)
        (let ((size (inexact size)))
          (new size 
               (make-rectangular (fl/ size 2.0) (fl/ size 2.0))
               (flsqrt 2.0)
               0.0
               (list)
               (random-seed)))))))

(define (in-bounds e position)           ;; Env Position -> Boolean
  ;; predicate: is |position| inside environment |e|?
  (let ((x (cfl-real-part position))
        (y (cfl-imag-part position)))
    (and (fl<=? 0.0 x) (fl<? x (env-size e))
         (fl<=? 0.0 y) (fl<? y (env-size e)))))
                                                                                            ;
(define (move e record-course)           ;; Env Boolean ->
  ;; update angle, position, and course of |e|
  (let ((prev-seed (random-seed)))       ;; decouple sequence of randoms from
    (random-seed (env-seed e))           ;; other uses of random
    (let ((max-rotation (fl/ (fl* 2.0 pi) 20.0)))
      (env-angle-set! e (fl+ (env-angle e)
                             (fl- max-rotation)
                             (random (fl* 2.0 max-rotation))))
      (let ((new-position
              (+ (env-position e)
                 (make-polar (env-speed e) (env-angle e)))))
        (cond
          [ (in-bounds e new-position)
              (env-position-set! e new-position)
              (env-seed-set! e (random-seed))
              (random-seed prev-seed)
              (when record-course
                (env-course-set! e (cons new-position (env-course e))))]
          [ else
              (env-angle-set! e (fl+ (env-angle e) (random pi) (fl/ pi 2.0)))
              (move e record-course)])))))

(define (_main args)                     ;; {KWarg} ->
#| Create environment and grid cell variant of spatial pooler using |args| |#
  (let* (
      [argument   (lambda (key default)
                    (let ([specified (assoc key args)])
                      (if specified (cdr specified) default)))]
      [train-time (argument 'train-time  1000)] ;; default for non-compile run
      [encoder-w  (argument 'encoder-w   75)]
      [encoder-n  (argument 'encoder-n   2500)]
      [encoder-r  (argument 'encoder-r   5)]
      [enc-num-samples (argument 'enc-num-samples  5)]
      [gc-num-samples  (argument 'gc-num-samples  15)]
      ;; Setup
      [env        (make-env 256)]
      [enc        (ce:make-ce encoder-w encoder-n encoder-r)]
      [gcm        (gc:make-gc `(
        [b1                              . 0.05]
        [b2                              . ,(fl/ 0.05 3.0)]
        [input-dimensions                . ,(list (ce:ce-n enc))]
        [column-dimensions               . (100)]
        [potential-pct                   . 0.95]
        [num-active-columns-per-inh-area . ,(int<- (* 0.2 100))]
        [syn-perm-inactive-dec           . ,(perm 0.001)]
        [syn-perm-active-inc             . ,(perm 0.005)]
        [syn-perm-connected              . ,(perm 0.250)]
        [stimulus-threshold              . 0]
        [boost-strength                  . 0.0]
        [global-inhibition               . #t]
        [potential-radius                . ,(ce:ce-n enc)]
        [wrap-around                     . #t]))])            
                                                                                            ;
  (let ((file "HTM-scheme/projects/KropffTreves2008_reproduction/experiment.data"))
                                                                                            ;
#;> (define (encode position)            ;; Position -> SDR
      (let ([nearest-position 
              (make-rectangular
                (fx+ encoder-r (exact (flround (cfl-real-part position))))
                (fx+ encoder-r (exact (flround (cfl-imag-part position)))))])
        (ce:encode enc nearest-position)))
                                                                                            ;
#;> (define (receptive-field bitx step   ;; Nat Nat (Coord Coord -> SDR) -> {(Vectorof Coord)}
               proc)
      ;; produce list of points for which bitx of result of proc for point is set
      (let ((size (exact (env-size env))))
        (let xloop ([x (fx1- size)] [rf (list (vector 0 0))])
          (if (fxnonnegative? x)
            (xloop (fx- x step)
              (let yloop ([y (fx1- size)] [rf rf])
                (if (fxnonnegative? y)
                    (let ([sdr (proc x y)])
                      (if (bitwise-bit-set? sdr bitx)
                        (yloop (fx- y step) (cons (vector x y) rf))
                        (yloop (fx- y step) rf)))
                    rf )))
            rf))))
                                                                                            ;
#;> (define (data-for-gcrf samp res)     ;; Nat Nat -> PlotData
      ;; produce gc receptive field plot data with resolution res
      [list "enc-rf"
        (receptive-field samp res (lambda (x y)
            (gc:reset gcm)
            (list->bitwise
              (gc:compute gcm (ce:encode enc (make-rectangular x y)) #f))))])
                                                                                            ;
#;> (define (data-for figure samp)       ;; Nat -> PlotData
      ;; produce plot data for one (sub)figure
      (let* (
          [data [list
            (case figure
              [("kt1")
                [list "course" (env-course env)]]
              [("kt2")
                [list "enc-rf"
                  (receptive-field samp 1 (lambda (x y)
                      (ce:encode enc (make-rectangular x y)))) ]]
              [("kt3" "kt4")
                  (data-for-gcrf samp 1)])]]
          [data (if samp  data
                  [cons [list "using" args] data])])
        [cons [list "figure" figure] data]))
                                                                                            ;
#| Training |#
    (gc:reset gcm)
    (let* (
        [gc-samples (vector->list
          (vector-sample (build-vector (sp:sp-num-columns gcm) id) gc-num-samples))]
        [enc-samples (take 5 gc-samples)]
        [untrained
          (map (lambda (s)
              (data-for "kt3" s))
            (take 5 gc-samples))] )
      (for-each display `("Training for " ,train-time " cycles ...\n"))
      (time
        (do-with-progress train-time (lambda (step)
            (when (zero? (fxmod step 50000))
              (with-output-to-file file (lambda ()
                  (write
                    (map (lambda (s)
                        [cons [list "figure" "kt4"]
                          [list (data-for-gcrf s 3)]])
                      (take 5 gc-samples))))
                'truncate)
              (when (fx=? step 50000)
                (ce:ce-radius-set! enc 6)
                (ce:ce-w-set! enc 100))
              (when (fx=? step 100000)
                (ce:ce-radius-set! enc 7)
                (ce:ce-w-set! enc 125))
              (when (fx=? step 150000)
                (ce:ce-radius-set! enc 8)
                (ce:ce-w-set! enc 150)))
            (move env (fx<? step 10000))
            (gc:compute gcm (encode (env-position env)) #t))))
                                                                                            ;
#| Save data for figures |#
      (with-output-to-file file (lambda ()
          (write
            (append
              [list (data-for "kt1" #f)]
              (map (lambda (s)
                  (data-for "kt2" s))
                enc-samples)
              untrained
              (map (lambda (s)
                  (data-for "kt4" s))
                gc-samples))))
        'truncate))  )))

(define (option name parameter)          ;; String [Number] -> KWarg
  ;; accept a few run options on command line
  (case name
    [("-tt")    `[train-time . ,(string->number parameter)] ]
    [("-ew")    `[encoder-w  . ,(string->number parameter)] ]
    [("-en")    `[encoder-n  . ,(string->number parameter)] ]
    [("-er")    `[encoder-r  . ,(string->number parameter)] ]
    [("-es")    `[enc-num-samples . ,(string->number parameter)] ]
    [("-gs")    `[gc-num-samples  . ,(string->number parameter)] ]
    ))
                                                                                            ;
(define (main)                           ;; {String} ->
  ;; eg $ scheme --program grid_cell_demo.wp -tt 10000
  (let parse ([cla (reverse (command-line-arguments))] [args (list)] [parameter "#t"])
    (cond
      [ (null? cla)  (_main args) ]
      [ (char=? #\- (string-ref (car cla) 0))
          (parse (cdr cla) (cons (option (car cla) parameter) args) "#t") ]
      [ else (parse (cdr cla) args (car cla)) ])))
        
(main)
        
  