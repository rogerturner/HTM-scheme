;; === HTM-scheme Layer4 algorithm  (C) 2021 Roger Turner. ===
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

*See htm_concept.ss for type and data structure description and code conventions*
  
Model layer 4 of a cortical column with inputs:
- tcs                thalamo-cortical ("feature")
- p6l4 (?)           ("location")
- ss4l4, ss4l23, p4  recurrent activity of this and nearby ccs  
Layer4 consists of ss4(L4), ss4(L23), and p4 populations, each modelled with the
Numenta apical_tiebreak_temporal_memory algorithm (apical input to p4 cells only).
The /main/ connections are shown below

               ^      ^       v
               |      |       |
               ^      ^       v
            +--ac-----pc-+-+--ai--+
      +---->bi |  p4  |  |*|      |
      |     +--|------|--|*|------+               ac = active cells
      |        ^      ^  |*|                      (output is combination
      |     +--ac-----pc-+*+------+                of p4 and ss4L23)
      +---->bi  ss4(L23) |*|      |
      |     +------------|*|------+               ai = apical input (from p2/3)
      |                  |*|
      |  +-----+         |*|                      bi = basal input
      |  |     ^         |*|
      |  |  +--ac--------|*|------+               pc = predicted cells
      |  +->bi  ss4(L4)  |*|      |
      |     +------------+-+------+
      |                   ^ \
      |                   |  \ minicolumn: any depolarized active cell inhibits others
      ^                   ^ 
 p6(L4) "location"    tcs "feature"
    
Connections are based on
[Izhikevich & Edelman 2008 "Large-scale model of mammalian thalamocortical systems"]
(10.1073/pnas.0712231105) figure 2, and figures 8 and 9 in Supporting Information.
Cells in each population receive input from multiple presynaptic populations
  |#

  #!chezscheme

(library (HTM-scheme HTM-scheme algorithms layer4)
                                                                                            ;
(export
make-l4
depolarize-cells
inhibit-columns
activate-cells
reset
get-active-cells
get-learning-cells
get-predicted-cells
get-inhibited-cols
l4-bursting-cols
  cols-from-cells
  reset-seq
  get-predicted-active-cells
  number-of-basal-segments
  print-statistics
  )
                                                                                            ;
(import
  (except (chezscheme) reset)
  (parameters)
  (HTM-scheme HTM-scheme algorithms htm_prelude)
  (HTM-scheme HTM-scheme algorithms htm_concept)
  (HTM-scheme HTM-scheme math       coordinates)
  (prefix (HTM-scheme HTM-scheme algorithms apical_tiebreak_temporal_memory) attm:))
                                                                                            ;
  (implicit-exports #f)

  ;; source layer coding:
  ;;    0  "location" (p6l4?)
  ;;    1  apical-input (p23)
  ;;    2  ss4l4
  ;;    3  ss4l23
  ;;    4  p4

;; === Layer record ===
                                                                                            ;
(define-record-type l4 (fields           ;; L4
    ss4l4                                ;; ATTM [ss4, axon -> L4]
    ss4l23                               ;; ATTM [ss4, axon -> L2/3]
    p4                                   ;; ATTM [p4, with apical dendrite]
    (mutable ss4l4-inhibited)            ;; {ColX}
    (mutable ss4l23-inhibited)
    (mutable p4-inhibited)
    (mutable bursting-cols))
(protocol #;(make-l4)
  ;; configure attm instances with connect? proc for 3 types of L4 cells
  (lambda (new)
    (lambda (
        ccx                              ;; index of this cortical column
        nmc                              ;; number of minicolumns/cc (127 if ncc>1)
        ss4l4-c/mc                       ;; cells/minicolumn for ss4l4 cells
        ss4l23-c/mc p4-c/mc 
        attm-overrides connect?)
;; (make-l4)
    (let ([dimensions `(
            [column-count      . ,nmc]
            [cortical-column   . ,ccx])])
      (new
        (attm:make-tm (append
            `([cells-per-column     . ,ss4l4-c/mc]
              #;[sample-size          . 8]
              [axon-radius2         . 1 #;500]
              [layer                . 2]
              [basal-connect?       . ,(lambda (source post-mcx)
                                         (connect? source post-mcx ccx 'ss4l4))])
            dimensions attm-overrides))
        (attm:make-tm (append
            `([cells-per-column     . ,ss4l23-c/mc]
              #;[sample-size          . 8]
              [axon-radius2         . 1 #;400]
              [permanence-decrement . ,(perm 0.025)]
              [layer                . 3]
              [basal-connect?       . ,(lambda (source post-mcx)
                                         (connect? source post-mcx ccx 'ss4l23))])
            dimensions attm-overrides))
        (attm:make-tm (append
            `([cells-per-column     . ,p4-c/mc]
              #;[sample-size          . 8]
              [axon-radius2         . 1 #;9]
              [permanence-decrement . ,(perm 0.03)]
              [layer                . 4]
              [basal-connect?       . ,(lambda (source post-mcx)
                                         (connect? source post-mcx ccx 'p4b))]
              [apical-connect?      . ,(lambda (source post-mcx)
                                         (connect? source post-mcx ccx 'p4a))])
            dimensions attm-overrides))
        (list) (list) (list) (list)))))))

;; === Layer algorithm ===
                                                                                            ;
(define (depolarize-cells l              ;; L4 {Source} {Source} {Source} Boolean ->
          location-input                 ;; (for this and neighbouring ccs)
          l4-activity                    ;; previous l4
          l4-apical-activity             ;; previous p23
          learn)
  ;; Depolarize each population based on current location input and previous layer activity
  (let ([basal-input  (sort-unique! (append location-input l4-activity '()))])
    (attm:depolarize-cells (l4-ss4l4  l) basal-input '() learn)
    (attm:depolarize-cells (l4-ss4l23 l) basal-input '() learn)
    (attm:depolarize-cells (l4-p4     l) basal-input l4-apical-activity learn)))
                                                                                            ;
(define (inhibit-columns l               ;; L4 {ColX} {Source} ->
          tcs location)
  ;; Calculate columns inhibited by early firing of predicted cells
  (define (inhibited-columns tm)         ;; TM -> {ColX}
    ;; Nearby cols projected to by predicted cells with TC input
    (if (fxzero? (attm:number-of-cells tm))  '()
      (let* (
          [correct-predicted-cells
            (attm:cells-in-cols tm (attm:get-predicted-cells tm) tcs)]
          [correct-predicted-cols
            (attm:cols-from-cells tm correct-predicted-cells)])
        (sort-unique! correct-predicted-cols))))
  (l4-ss4l4-inhibited-set!  l (inhibited-columns (l4-ss4l4  l)))
  (l4-ss4l23-inhibited-set! l (inhibited-columns (l4-ss4l23 l)))
  (l4-p4-inhibited-set!     l (inhibited-columns (l4-p4     l))))
                                                                                            ;
(define (activate-cells l                ;; L4 {Source} {Source} {Source} {Source} {Source} {Source} Boolean ->
          tcs                            ;; thalamo-cortical driving input "feature"
          location-input                 ;; context input (l6p4?)
          l4-activity                    ;; previous l4
          l4-learning                    ;; previous l4
          l4-apical-activity             ;; previous p23 + ss4l23 + p4
          l4-apical-learning             ;; ss4l23 + p4
          learn)
  ;; Activate and learn based on all active/learning cells connected to this cc and its inhibited cols
  (let (
      [basal-input       (sort-unique! (append location-input l4-activity '()))]
      [basal-growth      (sort-unique! (append location-input l4-learning '()))]
      [bursting-columns  (setdiff1d tcs
              (union1d (l4-ss4l4-inhibited l) (l4-ss4l23-inhibited l) (l4-p4-inhibited l)))])
    (attm:activate-cells (l4-ss4l4 l)
      tcs                        ;; proximal-input
      basal-input                ;; basal-reinforce-candidates
      '()                        ;; apical-reinforce-candidates
      basal-growth               ;; basal-growth-candidates
      '()                        ;; apical-growth-candidates
      learn 
      bursting-columns)
    (attm:activate-cells (l4-ss4l23 l)
      tcs                          ;; proximal-input
      basal-input                  ;; basal-reinforce-candidates
      '()                          ;; apical-reinforce-candidates
      basal-growth                 ;; basal-growth-candidates
      '()                          ;; apical-growth-candidates
      learn 
      bursting-columns)
    (attm:activate-cells (l4-p4 l)
      tcs                          ;; proximal-input
      basal-input                  ;; basal-reinforce-candidates
      l4-apical-activity           ;; apical-reinforce-candidates
      basal-growth                 ;; basal-growth-candidates
      l4-apical-learning           ;; apical-growth-candidates
      learn 
      bursting-columns)
    (l4-bursting-cols-set! l bursting-columns)))
                                                                                            ;
(define (reset l)                        ;; L4 ->
  (attm:reset (l4-ss4l4 l))
  (attm:reset (l4-ss4l23 l))
  (attm:reset (l4-p4 l)))
                                                                                            ;
(define (reset-seq l)                    ;; L4 ->
  (attm:reset (l4-ss4l4 l)))

;; === Accessors ===
                                                                                            ;
(define (get-inhibited-cols l pop)       ;; L4 L4Pop -> {ColX}
  ;; (use pop selector for compatibility)
  (case pop
    [ (ss4l4)  (l4-ss4l4-inhibited  l) ]
    [ (ss4l23) (l4-ss4l23-inhibited l) ]
    [ (p4)     (l4-p4-inhibited     l) ] ))
                                                                                            ;
(define (get f l pop)                    ;; (ATTM -> X) L4 L4Pop -> X
  ;; access f cells of sub-population pop of layer l
  (case pop
    [ (ss4l4)  (f (l4-ss4l4  l)) ]
    [ (ss4l23) (f (l4-ss4l23 l)) ]
    [ (p4)     (f (l4-p4     l)) ] ))
                                                                                            ;
(define (getx f l pop x)                 ;; (ATTM -> X) L4 L4Pop -> X
  ;; access f cells of sub-population pop of layer l
  (case pop
    [ (ss4l4)  (f (l4-ss4l4  l) x) ]
    [ (ss4l23) (f (l4-ss4l23 l) x) ]
    [ (p4)     (f (l4-p4     l) x) ] ))
                                                                                            ;
(define (get-active-cells l pop)         ;; L4 L4Pop -> {CellX}
  (get attm:get-active-cells l pop))
                                                                                            ;
(define (get-learning-cells l pop)       ;; L4 L4Pop -> {CellX}
  (get attm:get-learning-cells l pop))
                                                                                            ;
(define (get-predicted-cells l pop)      ;; L4 L4Pop -> {CellX}
  (get attm:get-predicted-cells l pop))
                                                                                            ;
(define (get-predicted-active-cells l p) ;; L4 L4Pop -> {CellX}
  (get attm:get-predicted-active-cells l p))
                                                                                            ;
(define (cols-from-cells l pop cellxs)   ;; L4 L4Pop {CellX} -> {ColX}
  (getx attm:cols-from-cells l pop cellxs))

;; === Statistics ===
                                                                                            ;
(define (number-of-connected-cells l p ) ;; L4 L4Pop -> Nat
  (get attm:number-of-connected-cells l p))
                                                                                            ;
(define (number-of-basal-synapses l pop) ;; L4 L4Pop -> Nat
  (get attm:number-of-basal-synapses l pop))
                                                                                            ;
(define (number-of-basal-segments l pop) ;; L4 L4Pop -> Nat
  (get attm:number-of-basal-segments l pop))
                                                                                            ;
(define (number-of-apical-synapses l p)  ;; L4 L4Pop -> Nat
  (get attm:number-of-apical-synapses l p))
                                                                                            ;
(define (number-of-apical-segments l p)  ;; L4 L4Pop -> Nat
  (get attm:number-of-apical-segments l p))
                                                                                            ;
(define (connection-lengths l p x)
  (getx attm:connection-lengths l p x))
                                                                                            ;
(define (print-statistics l4s)           ;; (Vector CCX->L4)
  (define (sum-l4 f pop)
    (vector-fold-left (lambda (sum layer)
        (+ sum (f layer pop)))
      0
      l4s))
  (let ([label "syns/segs/cells "])
    (define (print-pop pop n-syns n-segs)
      (let ([ncc (sum-l4 number-of-connected-cells pop)])
        (when (positive? ncc)
          (for-each display `( ,label
              ,(sum-l4 n-syns pop) "/"
              ,(sum-l4 n-segs pop) "/"
              ,ncc " " ,(symbol->string pop) "\n"))
          (set! label "                "))))
    (print-pop 'ss4l4 number-of-basal-synapses number-of-basal-segments)
    (print-pop 'ss4l23 number-of-basal-synapses number-of-basal-segments)
    (print-pop 'p4 number-of-basal-synapses number-of-basal-segments)
    (print-pop 'p4 number-of-apical-synapses number-of-apical-segments))
  (let ([label "pre-index (avg) "])
    (define (mean-l4 pop kind)
      (vector-fold-left (lambda (acc l4)
          (call-with-values
            (lambda () (connection-lengths l4 pop kind))
            (lambda (total count)
              (cons (+ (car acc) total) (+ (cdr acc) count)))))
        (cons 0 0)
        l4s))
    (define (print-pre-index pop kind)
      (for-each display `( ,label
          ,(let ([counts (mean-l4 pop kind)])
             (if (zero? (cdr counts))  0
                 (/ (inexact (quotient (* 10 (car counts)) (cdr counts))) 10)))
          " " ,(symbol->string pop) "\n"))
          (set! label "                "))
    (print-pre-index 'ss4l4  'basal)
    (print-pre-index 'ss4l23 'basal)
    (print-pre-index 'p4     'basal)
    (print-pre-index 'p4     'apical))
    
    (when (>= (vector-length l4s) 7) (print-max-axon-arbor l4s)))
                                                                                            ;
(define (print-max-axon-arbor l4s)
  ;; print minicolumn lattice showing post-synaptic cols for cell in cc0 with most synapses
  (define (segs->ccolx tm source segs)   ;; {Seg} -> {Nat}
    ;; map segs for a source to ccx+colx
    (let ([connected (attm:tm-connected-permanence tm)])
      (map (lambda (seg)
          (fx+ (fx* #x10000 (seg-ccx seg)) (attm:cellx->colx tm (seg-cellx seg))))
        (filter (lambda (seg)            ;; only segs that are connected
            (exists (lambda (synapse)
                (and (fx=? source (syn-source synapse))
                     (fx>=? (syn-perm synapse) connected)))
              (seg-synapses->list seg)))
          segs))))
  (let* (
      [ncc (min 37 (vector-length l4s))]
      [all-from-cc0                      ;; {[Source . {CColX}]}
        ;; list of pairs: [p4 cell in cc0 . list of cc+colx of projected to segments]
        (vector-fold-left                ;; accumulate over ccs
          (lambda (acc-per-cc l4 ccx)
            ;; accumulate each layer
            (define (from-p4-cc0 tm)    ;; TM -> {[Source . {CColX}]}
              ;; produce post colx from other ccs where source is p4 in cc0
              (let ([connected (attm:tm-connected-permanence tm)])
                (if (fxpositive? ccx)
                  (let-values ([(sources segss) (attm:get-axon-tree tm 'basal)])
                    (vector-fold-left    ;; accumulate over sources
                      (lambda (acc-per-source source segs)
                        ;; accumulate each source
                        (if (and (fx=? 0 (source-ccx source)) 
                                 (fx=? 2 (source-layer source))
                                 (pair? segs))
                          (cons (cons source (segs->ccolx tm source segs)) acc-per-source)
                          acc-per-source))
                      (list)
                      sources segss))
                  (list))))
            (append (from-p4-cc0 (l4-ss4l4 l4))
                    (from-p4-cc0 (l4-ss4l23 l4))
                    (from-p4-cc0 (l4-p4 l4))
                    acc-per-cc))
          (list)
          (vector-take ncc l4s) (list->vector (iota ncc))) ]
      [all-from-cc0                      ;; {[Source . {CColX}]}
        (sort (lambda (p1 p2) (fx<? (car p1) (car p2))) all-from-cc0)]
      [condensed-from-cc0                ;; {[Source . {CColX}]}
        ;; for each source, combine projected to lists for ccs
        (let each-run ([sources all-from-cc0] [out (list)])
          (cond
            [(null? sources) out ]
            [else
              (let ([this-source (caar sources)])
                (let each-source ([sources sources] [ccolxs (list)])
                  (cond
                    [(or (null? sources) (not (fx=? this-source (caar sources))))
                       (each-run sources (cons (cons this-source ccolxs) out)) ]
                    [else (each-source (cdr sources) (append (cdar sources) ccolxs)) ]))) ])) ]
      [max-from-cc0                      ;; [Source . {CColX}]
        ;; cell in cc0 with most connections to other ccs
        (fold-left
          (lambda (max-so-far source+posts)
            (if (fx>? (length source+posts) (length max-so-far))  source+posts
                max-so-far))
          '()
          condensed-from-cc0)]
      [max-from-cc0
        (if (pair? max-from-cc0)  max-from-cc0
            (list (make-source 0 2 0)))]
      [posts-in-cc0                      ;; {CColX}
        ;; post cols in cc0 for that cell
        (let ([source (car max-from-cc0)]
              [tm (l4-ss4l4 (vector-ref l4s 0))])
          (let-values ([(sources segss) (attm:get-axon-tree tm 'basal)])
            (let ([segs
                (do ([i 0 (fx1+ i)])
                    ((or (fx=? i (vector-length sources))
                         (fx=? source (vector-ref sources i)))
                      (vector-ref segss i)))])
              (segs->ccolx tm source segs))))]
      [max-from-cc0                      ;; [Source . {CColX}]
        (cons (car max-from-cc0)
          (append (cdr max-from-cc0) posts-in-cc0))])
    (define (print height width)          
      (let ([canvas (build-vector height (lambda _ (make-string width #\ )))])
        (define (row ccx mcx)
          (+ (fxdiv height 2) (r-coord-of-cc-centre ccx) (r-coord-of-minicol mcx)))
        (define (col ccx mcx)
          (+ (fxdiv width 2) (r-coord-of-cc-centre ccx) (r-coord-of-minicol mcx)
                (* 2 (q-coord-of-cc-centre ccx)) (* 2 (q-coord-of-minicol mcx))))
        (define (set-col! ccx mcx ch)
          (string-set!
            (vector-ref canvas (row ccx mcx))
            (col ccx mcx) ch))
        (do ([ccx 0 (+ ccx 1)]) ((= ccx ncc))
          (do ([mcx 0 (+ mcx 1)]) ((= mcx minicolumns/macrocolumn))
            (set-col! ccx mcx #\x22C5 )))
        (for-each (lambda (ccolx)
            (let* (
                [ccx  (fxdiv ccolx #x10000)]
                [mcx  (fxmod ccolx #x10000)]
                [ch   (string-ref (vector-ref canvas (row ccx mcx)) (col ccx mcx))])
              (set-col! ccx mcx
                (case ch
                  [(#\x22C5)  #\x2776 ]
                  [(#\x2776)  #\x2777 ]
                  [(#\x2777)  #\x2778 ]
                  [(#\x2778)  #\x2779 ]
                  [(#\x2779)  #\x277A ]
                  [else       #\x25CF ]))))
          (cdr max-from-cc0))
        (set-col! 0 (attm:cellx->colx (l4-ss4l4 (vector-ref l4s 0))
                      (source-cellx (car max-from-cc0))) #\x25C9)
        (for-each display `( ,(length max-from-cc0) " post-synaptic minicolumns\n"))
        (vector-for-each (lambda (s)
            (put-string (current-output-port) s) (newline))
          canvas)))
    (let ([radius (case minicolumns/macrocolumn
                    [(127)  6]
                    [(169)  7]
                    [(217)  8]
                    [(271)  9]
                    [(331) 10]
                    [(397) 11])])
      (print (case ncc
              [(7)  (fx+ 3 (fx* radius 6))]
              [(19) (fx+ 5 (fx* radius 10))]
              [(37) 91])
             (case ncc
              [(7)  (fx- (fx* radius 12) 6)]
              [(19) (fx- (fx* radius 20) 15)]
              [(37) 145])))))

)