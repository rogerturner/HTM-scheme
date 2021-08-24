#!chezscheme

;; === HTM-scheme Untangling Sequences experiment  (C) 2019-2021 Roger Turner. ===
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Based on code from Numenta Platform for Intelligent Computing (NuPIC) ;;
  ;; which is Copyright (C) 2017, Numenta, Inc.                            ;;
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

  ;; Partial replication and extension of experiments reported in Numenta paper
  ;; Ahmad & Hawkins 2017 "Untangling Sequences: Behavior vs. External Causes"
  ;; doi: 10.1101/190678 
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.
  ;;
  ;; "Remember that all models are wrong; the practical question is how wrong
  ;;  do they have to be to not be useful" [George Box]

(import
  (chezscheme)
  (untangling_sequences-lib))

(define exp4a '(
    [figure         .  f4a]
    [num-sequences  .    5]
    [num-objects    .    0]
    [num-features   .   10]
    [num-locations  .  100]
    [num-repetitions .   20]
    [train-keys     . (l4ac TMac TXac)]
    [test-keys      . (l4lnp l4pa TMlnp TMpa TXlnp TXpa)]))
                                                                                            ;
(define exp5a '(
    [figure         .  f5a]
    [num-objects    .   50]
    [num-sequences  .    0]
    [num-features   .  100]
    [num-locations  .  100]
    [train-keys     . (l4ac TMac TXac)]
    [test-keys      . (l4lnp l4pa TMlnp TMpa TXlnp TXpa)]))
                                                                                            ;
(define exp6 '(                          ;; params from combined_sequences.py runExperiment6
    [figure          .  f6]
    [num-sequences   .  50]
    [num-objects     .  50]
    [num-features    .  50]
    [num-locations   . 100]
    [num-repetitions .  32]
    [interleave-training   . #t]
    [random-seq-location   . #t]
    [location-per-sequence . #f]
    [train-keys      . (l4ac TMac TXac #;l3ac)]
    [test-keys       . (l4pa TMpa TXpa #;l3pa) #;(l4lpa TMlpa TXlpa)]))
                                                                                            ;
(define run                              ;; String [ {KWarg} ] ->
  ;; can be used in repl eg (run "4a" '( [column-count . 150] ) )
  (case-lambda 
    [ (figure) (run figure '()) ]
    [ (figure options)
        (let ([start (statistics)])
          (collect) (collect) (collect (collect-maximum-generation) 'static)
          (experiment
            (case figure
              [("4a")    exp4a ]
              [("5a")    exp5a ]
              [("6")     exp6  ])
            options)
          (collect) (collect) (collect) (collect)
          (sstats-print (sstats-difference (statistics) start))) ] ))
                                                                                            ;
(define (option name parameter)          ;; String [Number | "#f" | "f" | X] -> KWarg
  ;; accept a few run options on command line
  (define (number) (string->number parameter))
  (define (boolean)
    (case parameter
      [("f" "#f") #f]
      [else #t]))
  (case name
    [("-nmc")   (cons 'column-count          (number))]
    [("-nib")   (cons 'num-input-bits        (number))]
    [("-l4")    (cons 'ss4l4-cells/mcol      (number))]
    [("-l23")   (cons 'ss4l23-cells/mcol     (number))]
    [("-p4")    (cons 'p4-cells/mcol         (number))]
    [("-in")    (cons 'intersperse-noise     (number))]
    [("-it")    (cons 'interleave-training   (boolean))]
    [("-lps")   (cons 'location-per-sequence (boolean))]
    [("-ncc")   (cons 'num-cortical-columns  (number))]
    [("-nf")    (cons 'num-features          (number))]
    [("-nl")    (cons 'num-locations         (number))]
    [("-no")    (cons 'num-objects           (number))]
    [("-nr")    (cons 'num-repetitions       (number))]
    [("-ns")    (cons 'num-sequences         (number))]
    [("-ol")    (cons 'online-learning       (boolean))]
    [("-rsl")   (cons 'random-seq-location   (boolean))]
    [("-ubc")   (cons 'use-bursting-columns  (boolean))]
    [("-ss")    (cons 'superimpose-sequence  (boolean))]))
                                                                                            ;
(define (main command-line)              ;; {String} ->
  ;; eg $ scheme --program untangling_sequences.wp 6 -cc 150 -sc f
  (let process ([args (reverse command-line)] [options (list)] [parameter "t"])
    (cond
      [(or (null? args)
           (string=? "" (car args))) ]   ;; (not run as --program)
      [(null? (cdr args))
        (run (car args) options)]       ;; last (first) element is experiment
      [(char=? #\- (string-ref (car args) 0))
        (process (cdr args) (cons (option (car args) parameter) options) "t")]
      [else
        (process (cdr args) options (car args))])))
                                                                                            ;
(let ([seed 1629305071 #;(time-second (current-time))])
  (random-seed seed)  ;; 1629179167 1629179375 1629295454 1629295597 1629305071
  (display seed) (display "  "))
        
(main (command-line-arguments))

#;(profile-dump-data "untangling_sequences.data") 

