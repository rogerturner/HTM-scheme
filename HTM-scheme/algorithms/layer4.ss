#!chezscheme

;; === HTM-scheme Layer4 algorithm Copyright 2019 Roger Turner. ===
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

  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.
  #|
  
Layer4 consists of ss4(L4), ss4(L23), and p4 populations, each modelled with the
Numenta apical_tiebreak_temporal_memory algorithm (apical input to p4 cells only).
The main connections are shown below (see figure 2 in Izhikevich & Edelman 2008 
"Large-scale model of mammalian thalamocortical systems" 10.1073/pnas.0712231105)

               ^      ^       v
               |      |       |
               ^      ^       v
            +--ac-----pc-+-+--ai--+
      +---->bi |  p4  |  |*|      |
      |     +--|------|--|*|------+               ac = active cells
      |        ^      ^  |*|                      (output is combination
      |     +--ac-----pc-+*+------+                of p4 and ss4L23)
      +---->bi  ss4(L23) |*|      |
      |     +------------|*|------+               ai = apical input
      |                  |*|
      +--------+         |*|                      bi = basal input
      |        ^         |*|                      (combination of location
      |     +--ac--------|*|------+                and ss4L4's active cells
      +---->bi  ss4(L4)  |*|      |                on previous step)
      |     +------------+-+------+
      |                   ^ \                     pc = predicted cells
      |                   |  \
      ^                   ^   \
    location        feature   minicolumn: any depolarized active cell inhibits others

Cells in each population receive input from multiple presynaptic populations (see
figures 8 and 9 in the Supporting Information for Izhikevich & Edelman 2008).

    |#

(library (HTM-scheme HTM-scheme algorithms layer4)
                                                                                            ;
(export
  l4-pop
  l4-active-cells
  l4-predicted-cells
  l4-pre-stride
  make-l4
  compute
  reset
  reset-seq
  get-active-cells
  get-predicted-cells
  get-predicted-active-cells
  number-of-basal-synapses
  number-of-basal-segments
  number-of-apical-synapses
  number-of-apical-segments)
                                                                                            ;
(import
  (except (chezscheme) reset)
          (HTM-scheme HTM-scheme algorithms htm_prelude)
          (HTM-scheme HTM-scheme algorithms htm_concept)
  (prefix (HTM-scheme HTM-scheme algorithms apical_tiebreak_temporal_memory) attm:))
                                                                                            ;
(define-enumeration l4-populations (ss4l4 ss4l23 p4) l4-pop)
                                                                                            ;
(define-record-type l4 (fields           ;; L4
  ss4l4                                  ;; ATTM ss4, axon -> L4
  ss4l23                                 ;; ATTM ss4, axon -> L2/3
  p4                                     ;; ATTM p4
  (mutable active-cells)
  (mutable predicted-cells)
  pre-stride)                            ;; Nat
  (protocol
    (lambda (new)
      (lambda (cc ncl4pop nib bis ais attm-overrides)
        (let* (
            (pre-stride (fxmax bis ais (fx* cc ncl4pop)))
            (dimensions `(
              [column-count      . ,cc]
              [cells-per-column  . ,ncl4pop]
              [sample-size       . ,(fx* 5 nib)]
              [basal-input-size  . ,(fx* 5 pre-stride)]
              [apical-input-size . ,(fx* 5 pre-stride)])))
          (new
            (attm:make-tm (append 
              '([has-recurrent-connections . #t]) dimensions attm-overrides))
            (attm:make-tm (append dimensions attm-overrides))
            (attm:make-tm (append dimensions attm-overrides))
            '()
            '()
            pre-stride))))))
                                                                                            ;
(define (compute l feature location      ;; L4 SDR SDR SDR Boolean -> SDR SDR SDR
          apical-input learn)
  ;; step layer l and return active non-bursting and predicted L2 projecting cells and bursting columns
  (define (offset inputs source)         ;; SDR Nat -> SDR
    ;; produce inputs with elements offset by (* source pre-stride)
    (if (zero? source)  inputs
      (let ((stride (fx* source (l4-pre-stride l))))
        (map (lambda (e) (fx+ e stride))
          inputs))))
  (define (thin inputs n)                ;; SDR Nat -> SDR
    ;; produce selection from inputs, where each element has 1/n chance of inclusion
    (if (fx=? n 1)  inputs
      (let ((v (list->vector inputs)))
        (sort! fx<? (vector->list (vector-sample v (int<- (/ (vector-length v) n))))))))
  (let* ( (ss4l4->        (attm:get-active-cells (l4-ss4l4 l)))
          (ss4l23->       (attm:get-active-cells (l4-ss4l23 l)))
          (p4->           (attm:get-active-cells (l4-p4 l)))
          (->l4-basal     (append (offset location              0)
                                  (offset (thin apical-input 6) 1)
                                  (offset (thin ss4l4->      1) 2)
                                  (offset (thin ss4l23->     3) 3)
                                  (offset (thin p4->         3) 4)))
          (->ss4l4-basal  (append (offset location               0)
                                  (offset (thin apical-input 12) 1)
                                  (offset (thin ss4l4->       1) 2)
                                  (offset (thin ss4l23->      3) 3)
                                  (offset (thin p4->          3) 4)))
          (->p4-apical    (append (offset apical-input      1)
                                  (offset (thin ss4l23-> 8) 3)
                                  (offset (thin p4->     8) 4))))
    (attm:depolarize-cells (l4-ss4l4  l) ->ss4l4-basal '()         learn)
    (attm:depolarize-cells (l4-ss4l23 l) ->l4-basal    '()         learn)
    (attm:depolarize-cells (l4-p4     l) ->l4-basal    ->p4-apical learn)
    (let* (
        (ss4l4-predicted-cells  (attm:get-predicted-cells (l4-ss4l4 l)))
        (ss4l4-predicted-cols   (attm:cols-from-cells (l4-ss4l4 l) ss4l4-predicted-cells))
        (ss4l23-predicted-cells (attm:get-predicted-cells (l4-ss4l23 l)))
        (ss4l23-predicted-cols  (attm:cols-from-cells (l4-ss4l23 l) ss4l23-predicted-cells))
        (p4-predicted-cells     (attm:get-predicted-cells (l4-p4 l)))
        (p4-predicted-cols      (attm:cols-from-cells (l4-p4 l) p4-predicted-cells))
        (bursting-columns       (setdiff1d feature
                                    (union1d ss4l4-predicted-cols 
                                      (union1d ss4l23-predicted-cols p4-predicted-cols))))
        (active-non-bursters    (lambda (pop)
                                  (setdiff1d 
                                    (attm:get-active-cells (pop l))
                                    (attm:get-all-cells-in-columns (pop l) bursting-columns)))))
      (attm:activate-cells (l4-ss4l4 l)
        feature                          ;; feedforward-input
        ->ss4l4-basal                    ;; basal-reinforce-candidates
        '()                              ;; apical-reinforce-candidates
        (append (thin location 5)        ;; basal-growth-candidates
                (offset (thin (attm:get-winner-cells (l4-ss4l4 l))  2) 2)
                (offset (thin (attm:get-winner-cells (l4-ss4l23 l)) 8) 3)
                (offset (thin (attm:get-winner-cells (l4-p4 l))     7) 4))
        '()                              ;; apical-growth-candidates
        learn 
        bursting-columns)
      (attm:activate-cells (l4-ss4l23 l)
        feature                          ;; feedforward-input
        ->l4-basal                       ;; basal-reinforce-candidates
        '()                              ;; apical-reinforce-candidates
        (append (thin location 1)        ;; basal-growth-candidates
                (offset (thin (attm:get-winner-cells (l4-ss4l4 l))  2) 2)
                (offset (thin (attm:get-winner-cells (l4-ss4l23 l)) 8) 3)
                (offset (thin (attm:get-winner-cells (l4-p4 l))     7) 4))
      '()                              ;; apical-growth-candidates
        learn 
        bursting-columns)
      (attm:activate-cells (l4-p4 l)
        feature                          ;; feedforward-input
        ->l4-basal                       ;; basal-reinforce-candidates
        ->p4-apical                      ;; apical-reinforce-candidates
        (append (thin location 1)        ;; basal-growth-candidates
                (offset (thin (attm:get-winner-cells (l4-ss4l4 l))  2) 2)
                (offset (thin (attm:get-winner-cells (l4-ss4l23 l)) 8) 3)
                (offset (thin (attm:get-winner-cells (l4-p4 l))     7) 4))
        ->p4-apical                      ;; apical-growth-candidates
        learn 
        bursting-columns)
      (l4-active-cells-set! l
        (append (offset (active-non-bursters l4-ss4l23) 0) (offset (active-non-bursters l4-p4) 1)))
      (l4-predicted-cells-set! l
        (append (offset ss4l23-predicted-cells 0) (offset p4-predicted-cells 1)))
      #;(values
        (append (offset (active-non-bursters l4-ss4l23) 0) (offset (active-non-bursters l4-p4) 1))
        (append (offset ss4l23-predicted-cells 0) (offset p4-predicted-cells 1))
        bursting-columns))))
                                                                                            ;
(define (reset l)                        ;; L4 ->
  (attm:reset (l4-ss4l4 l))
  (attm:reset (l4-ss4l23 l))
  (attm:reset (l4-p4 l)))
                                                                                            ;
(define (reset-seq l)                    ;; L4->
  (attm:reset (l4-ss4l4 l)))
                                                                                            ;
(define (get f l pop)                    ;; (ATTM -> X) L4 L4Pop -> X
  ;; access f cells of sub-population pop of layer l
  (cond
    [ (enum-set-member? 'ss4l4 pop)  (f (l4-ss4l4  l)) ]
    [ (enum-set-member? 'ss4l23 pop) (f (l4-ss4l23 l)) ]
    [ (enum-set-member? 'p4 pop)     (f (l4-p4     l)) ] ))
                                                                                            ;
(define (get-active-cells l pop)         ;; L4 L4Pop -> {CellX}
  (get attm:get-active-cells l pop))
                                                                                            ;
(define (get-predicted-cells l pop)      ;; L4 L4Pop -> {CellX}
  (get attm:get-predicted-cells l pop))
                                                                                            ;
(define (get-predicted-active-cells l pop)
  (get attm:get-predicted-active-cells l pop))
                                                                                            ;
(define (number-of-basal-synapses l pop) ;; L4 L4Pop -> Nat
  (get attm:number-of-basal-synapses l pop))
                                                                                            ;
(define (number-of-basal-segments l pop)
  (get attm:number-of-basal-segments l pop))
                                                                                            ;
(define (number-of-apical-synapses l pop)
  (get attm:number-of-apical-synapses l pop))
                                                                                            ;
(define (number-of-apical-segments l pop)
  (get attm:number-of-apical-segments l pop))
                                                                                            ;
)