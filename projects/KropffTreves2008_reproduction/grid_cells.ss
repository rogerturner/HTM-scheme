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
  
Extends Spatial Pooler: see
(https://github.com/ctrl-z-9000-times/KropffTreves2008_reproduction)
  |#

(library (HTM-scheme projects KropffTreves2008_reproduction grid_cells)

(export
  make-gc
  reset
  compute)
  
(import
  (except (chezscheme) make-list reset)
          (HTM-scheme HTM-scheme algorithms htm_prelude)
          (HTM-scheme HTM-scheme algorithms htm_concept)
  (prefix (HTM-scheme HTM-scheme algorithms spatial_pooler) sp:))

(define-record-type gc
  (parent sp:sp)
  (fields
    b1
    b2
    (mutable r-act)
    (mutable r-inact)
    (mutable prev-active-columns)
    (mutable prev-input-vector))
  (protocol
    (lambda (pargs->new)                 ;; (listof KWarg) -> TM
      (lambda (kwargs)
        (apply (pargs->new kwargs)
               (key-word-args kwargs gc-defaults))))))
                                                                                            ;
(define gc-defaults `(
    [b1                  . 0.05]
    [b2                  . ,(/ 0.05 3)]
    [r-act               . #()]
    [r-inact             . #()]
    [prev-active-columns . () ]
    [prev-input-vector   . 0  ]))
    
(define (reset gc)
  ;; 
  (gc-r-act-set!   gc (make-vector (sp:sp-num-columns gc) 0.0))
  (gc-r-inact-set! gc (make-vector (sp:sp-num-columns gc) 0.0))
  (gc-prev-active-columns-set! gc (list))
  (gc-prev-input-vector-set!   gc 0 ))
                                                                                            ;
(define (compute gc input-vector learn)  ;; GC InputVec Boolean -> (listof ColumnX)
  ;; override calculate-overlap and adapt-synapses in sp:compute
  (sp:compute gc input-vector learn 
    (lambda (calculate-overlap)
      (let* (
          [overlaps         (calculate-overlap gc input-vector)]
          [overlaps         (vector-map (lambda (ov cs)
                                (let ([sum-synapses (bitwise-bit-count cs)])
                                  (if (fxzero? sum-synapses)  0.0 
                                      (fl/ ov (fixnum->flonum sum-synapses)))))
                              overlaps (sp:sp-connected-synapses gc))])
        (if learn (_fatigue gc overlaps)
                  overlaps)))
    (lambda (active-columns)
      (adapt-synapses gc input-vector active-columns))))
                                                                                            ;
(define (_fatigue gc overlaps)           ;; (ColVecOf Number) (ColVecOf OverlapX) -> (ColVecOf OverlapX)
  ;;
  (gc-r-act-set! gc 
    (vector-map (lambda (ra ri ov)
        (fl+ ra (fl* (gc-b1 gc) (fl- ov ri ra))))
      (gc-r-act gc) (gc-r-inact gc) overlaps))
  (gc-r-inact-set! gc 
    (vector-map (lambda (ri ov)
        (fl+ ri (fl* (gc-b2 gc) (fl- ov ri))))
      (gc-r-inact gc) overlaps))
  (gc-r-act gc))
                                                                                            ;
(define (adapt-synapses                  ;; SP InputVec (listof ColumnX) ->
          gc input-vector active-columns)
  ;; update permanences in segments of active columns (+ if synapse's input on, - if not)
  (let ((syn-perm-trim-threshold (sp:sp-syn-perm-trim-threshold gc))
        (syn-perm-active-inc     (sp:sp-syn-perm-active-inc gc))
        (syn-perm-inactive-dec   (sp:sp-syn-perm-inactive-dec gc))
        (syn-perm-connected      (sp:sp-syn-perm-connected gc)))
    (threaded-vector-for-each (lambda (cx)
        (let ([prev-active (memv cx (gc-prev-active-columns gc))]
              [segment (vector-ref (sp:sp-potential-pools gc) cx)])
          (if (zero? (sp:sp-stimulus-threshold gc))
            (do ([i (fx1- (synapses-length segment)) (fx1- i)])
                ((fxnegative? i))
              (let* ( [synapse (synapses-ref segment i)]
                      [inpx    (syn-prex synapse)]
                      [inp-bit (bitwise-bit-set? input-vector inpx)] )
                (unless (and prev-active ;; skip repeated permanence change
                          (eq? inp-bit (bitwise-bit-set? (gc-prev-input-vector gc) inpx)))
                  (let ([perm (syn-perm synapse)])
                    (cond
                      [(fx<? min-perm perm max-perm)
                        (synapses-set! segment i
                          (if inp-bit
                            (sp:increase-perm synapse syn-perm-active-inc)
                            (sp:decrease-perm synapse syn-perm-inactive-dec syn-perm-trim-threshold)))]
                      [(fx=? perm min-perm)
                        (when inp-bit
                          (synapses-set! segment i (make-syn inpx syn-perm-trim-threshold)))]
                      [else
                        (unless inp-bit
                          (synapses-set! segment i (make-syn inpx (fx- max-perm syn-perm-inactive-dec))))])))))
            (let loop ((i (fx1- (synapses-length segment))) (num-connected 0))
              (cond 
                [(fxnegative? i)
                  (when (fx<? num-connected (sp:sp-stimulus-threshold gc))
                    (sp:raise-permanence-to-threshold gc segment)) ]
                [else
                  (let* ( [synapse (synapses-ref segment i)]
                          [inpx    (syn-prex synapse)]
                          [inp-bit (bitwise-bit-set? input-vector inpx)] )
                    (if (and prev-active ;; skip repeated permanence change
                          (eq? inp-bit (bitwise-bit-set? (gc-prev-input-vector gc) inpx)))
                      (loop (fx1- i) (if (fx>=? (syn-perm synapse) syn-perm-connected)
                                          (fx1+ num-connected)
                                          num-connected))
                      (let ([synapse (if inp-bit
                                         (sp:increase-perm synapse syn-perm-active-inc)
                                         (sp:decrease-perm synapse syn-perm-inactive-dec syn-perm-trim-threshold))])
                        (synapses-set! segment i synapse)
                        (loop (fx1- i) (if (fx>=? (syn-perm synapse) syn-perm-connected)
                                            (fx1+ num-connected)
                                            num-connected)))))])))
          )
          (vector-set! (sp:sp-connected-synapses gc) cx (sp:connected-inputs gc cx)))
      (list->vector active-columns) #;(vector-sample (list->vector active-columns) 7)))
  (gc-prev-active-columns-set! gc active-columns)
  (gc-prev-input-vector-set!   gc input-vector))
                                                                                            ;
)
