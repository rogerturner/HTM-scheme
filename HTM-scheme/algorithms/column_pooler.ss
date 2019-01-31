#!chezscheme

;; === HTM-scheme Column Pooler Copyright 2019 Roger Turner. ===
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Based on code from Numenta Platform for Intelligent Computing (NuPIC) ;;
  ;; which is Copyright (C) 2016, Numenta, Inc.                            ;;
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

  ;; Translated from numenta htmresearch/.../column_pooler.py --
  ;; see comments there for descriptions of functions and parameters.
  ;; "This class constitutes a temporary implementation for a cross-column pooler.
  ;;  The implementation goal of this class is to prove basic properties before
  ;;  creating a cleaner implementation." [Numenta description]
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.

(library (HTM-scheme HTM-scheme algorithms column_pooler)
                                                                                            ;
(export
  make-cp
  compute
  reset
  (rename
    (cp-input-width  number-of-inputs)
    (cp-cell-count   number-of-cells)
    (cp-sdr-size     get-sdr-size)
    (cp-active-cells get-active-cells))
  (rename
    (num-connected-proximal-synapses test:num-connected-proximal-synapses)
    (number-of-proximal-synapses     test:number-of-proximal-synapses)
    (cp-proximal-permanences         test:cp-proximal-permanences)))
                                                                                            ;
(import 
  (except (chezscheme) add1 make-list random reset)
  (HTM-scheme HTM-scheme algorithms htm_prelude)
  (HTM-scheme HTM-scheme algorithms htm_concept))

;; === Layer record ===
                                                                                            ;
(define-record-type cp                   ;; CP
  (fields
    input-width
    lateral-input-widths
    cell-count
    sdr-size
    online-learning
    (mutable max-sdr-size)
    (mutable min-sdr-size)
    syn-perm-proximal-inc                ;; Proximal
    syn-perm-proximal-dec
    initial-proximal-permanence
    sample-size-proximal
    min-threshold-proximal
    connected-permanence-proximal
    predicted-inhibition-threshold
    syn-perm-distal-inc                  ;; Distal
    syn-perm-distal-dec
    initial-distal-permanence
    sample-size-distal
    activation-threshold-distal
    connected-permanence-distal
    inertia-factor
    seed
    (mutable active-cells)
    (mutable proximal-permanences)       ;; CellVecOf Segment
    (mutable internal-distal-permanences);; CellVecOf Segment
    (mutable distal-permanences)         ;; {CellVecOf Segment}
    use-inertia)
  (protocol
    (lambda (new)
      (lambda (kwargs)                   ;; (listof KWarg) -> CP
        (let ((cp (apply new (key-word-args kwargs cp-defaults))))
          (when (zero? (cp-max-sdr-size cp))
            (cp-max-sdr-size-set! cp (cp-sdr-size cp)))
          (when (zero? (cp-min-sdr-size cp))
            (cp-min-sdr-size-set! cp (cp-sdr-size cp)))
          (random-seed! (cp-seed cp))
          (cp-proximal-permanences-set! cp
            (make-vector (cp-cell-count cp) (fxvector)))
          (cp-internal-distal-permanences-set! cp
            (make-vector (cp-cell-count cp) (fxvector)))
          ;; (distal permanences are allocated when lateral inputs known)
          cp)))))
                                                                                           ;
(define cp-defaults `(                   ;; (listof KWarg)
    [input-width                    . 0]
    [lateral-input-widths           . ()]
    [cell-count                     . 4096]
    [sdr-size                       . 40]
    [online-learning                . #f]
    [max-sdr-size                   . 0]
    [min-sdr-size                   . 0]
    [syn-perm-proximal-inc          . ,(perm 0.1)]
    [syn-perm-proximal-dec          . ,(perm 0.001)]
    [initial-proximal-permanence    . ,(perm 0.6)]
    [sample-size-proximal           . 20]
    [min-threshold-proximal         . 10]
    [connected-permanence-proximal  . ,(perm 0.50)]
    [predicted-inhibition-threshold . 20]
    [syn-perm-distal-inc            . ,(perm 0.1)]
    [syn-perm-distal-dec            . ,(perm 0.001)]
    [initial-distal-permanence      . ,(perm 0.6)]
    [sample-size-distal             . 20]
    [activation-threshold-distal    . 13]
    [connected-permanence-distal    . ,(perm 0.50)]
    [inertia-factor                 . 1.]
    [seed                           . 42]
    [active-cells                   . ()]
    [proximal-permanences           . #()]
    [internal-distal-permanences    . #()]
    [distal-permanences             . ()]
    [use-inertia                    . #t]))

;; === Synapses and Permanences ===
;; ('Permanences' in column_pooler.py correspond to Synapses in HTM-scheme)
                                                                                            ;
(define (prune-synapses synapses omits)  ;; Synapses (listof Nat) -> Synapses
  ;; omit from synapses elements indexed by omits (which is sorted)
  (let* ( (reslen (fx- (synapses-length synapses) (length omits)))
          (result (make-synapses reslen)))
    (let loop ((rx 0) (sx 0) (omits omits))
      (cond [ (fx=? rx reslen) result ]
            [ (if (null? omits) #f
                (fx=? sx (car omits))) (loop rx (add1 sx) (cdr omits)) ]
            [ else
                (synapses-set! result rx (synapses-ref synapses sx))
                (loop (add1 rx) (add1 sx) omits) ]))))
                                                                                            ;
(define (adapt-synapses! permanences     ;; Permanences {CellX} {InputX} Perm Perm ->
          active-cells active-input permanence-increment permanence-decrement)
  ;; update synapses: strengthen those connected to input, weaken others,
  ;; remove synapse on zero permanence
  (for-each
    (lambda (cellx)
      (vector-set! permanences cellx 
        (let ((synapses (vector-ref permanences cellx))
              (active-input (list->vector active-input)))
          (let build-s-t-d [ (sx (fx- (synapses-length synapses) 1))
                             (synapses-to-destroy '()) ]
            (if (negative? sx)
                (cond [ (null? synapses-to-destroy) synapses ]
                      [ (fx=? (length synapses-to-destroy) (synapses-length synapses)) (make-synapses 0) ]
                      [ else (prune-synapses synapses synapses-to-destroy) ])
                (let* ( (synapse    (synapses-ref synapses sx))
                        (prex       (syn-prex synapse))
                        (permanence (if (fxsearch active-input prex)
                                      (clip-max (fx+ (syn-perm synapse) permanence-increment))
                                      (clip-min (fx- (syn-perm synapse) permanence-decrement)))))
                  (if (zero? permanence)
                      ;; build synapses-to-destroy indices as sorted list
                      (build-s-t-d (fx- sx 1) (cons sx synapses-to-destroy))
                      (begin
                        (synapses-set! synapses sx (make-syn prex permanence))
                        (build-s-t-d (fx- sx 1) synapses-to-destroy)))))))))
    active-cells))
                                                                                            ;
(define (add-synapses! permanences       ;; Permanences {CellX} {InputX} Perm ->
          active-cells active-input initial-permanence)
  ;; create new synapses for all input bits ('set-zeros-on-outer')
  (for-each
    (lambda (cellx)
      (let ((synapses (vector-ref permanences cellx)))
        (let add-synapses [ (all-synapses (fxvector->list synapses)) 
                            (inputs active-input) ]
          (if (null? inputs)
            (vector-set! permanences cellx (list->fxvector all-synapses))
            (let* ( (syn-low  (make-syn (car inputs) min-perm))
                    (syn-high (fx+ syn-low max-perm)))
              (if (synapses-search synapses syn-low syn-high)
                (add-synapses all-synapses (cdr inputs))
                (add-synapses 
                  (cons (make-syn (car inputs) initial-permanence) all-synapses)
                  (cdr inputs))))))
        (vector-set! permanences cellx
          (list->fxvector (list-sort fx<? (fxvector->list (vector-ref permanences cellx)))))))
    active-cells))
                                                                                            ;
(define (grow-synapses! permanences      ;; Permanences {CellX} {InputX} {Nat} Perm ->
          active-cells growth-candidates max-news initial-permanence)
  ;; create new synapses for some growth-candidates ('set-random-zeros-on-outer')
  (for-each
    (lambda (cellx max-new)
      (when (positive? max-new)
        (let ((synapses (vector-ref permanences cellx)))
          (let choose-new-synapses [ (new-synapses (list)) 
                                     (inputs growth-candidates) ]
            (if (null? inputs)
              (vector-set! permanences cellx
                (list->fxvector 
                  (append
                    (fxvector->list synapses)
                    (if (fx<? max-new (length new-synapses))
                        (vector->list (vector-sample (list->vector new-synapses) max-new))
                        new-synapses))))
              (let* ( (syn-low  (make-syn (car inputs) min-perm))
                      (syn-high (fx+ syn-low max-perm)))
                (if (synapses-search synapses syn-low syn-high)
                  (choose-new-synapses new-synapses (cdr inputs))
                  (choose-new-synapses 
                    (cons (make-syn (car inputs) initial-permanence) new-synapses)
                    (cdr inputs)))))))
        (vector-set! permanences cellx
          (list->fxvector (list-sort fx<? (fxvector->list (vector-ref permanences cellx)))))))
    active-cells max-news))
                                                                                            ;
(define (max-new-by-cell permanences     ;; Permanences {CellX} {InputX} Nat -> {Nat}
          active-cells active-input sample-size)
  ;; produce (sample-size - 'n-non-zeros-per-row-on-cols') for each cell
  (map
    (lambda (cellx)
      (let* [ (synapses (vector-ref permanences cellx))
              (num-synapses (synapses-length synapses)) ]
        (do [ (sx 0 (add1 sx))
              (count 0 (if (memv (syn-prex (synapses-ref synapses sx)) active-input)
                           (add1 count)
                           count)) ]
            ((fx=? sx num-synapses) (fx- sample-size count)))))
    active-cells))
                                                                                            ;
(define (right-vec-sum-at-nzgte-thresh   ;; Permanences {InputX} Perm -> (CellVecOf Nat)
          permanences inputxs threshold)
  ;; produce counts of overlaps between synapses above threshold and inputs
  (let* ( (inputxs   (list->vector inputxs))
          (syn-lows  (vector-map
                       (lambda (inputx) (make-syn inputx min-perm))
                       inputxs))
          (syn-highs (vector-map
                       (lambda (syn-low) (fx+ syn-low max-perm))
                       syn-lows)))
    (vector-map
      (lambda (synapses)
        (vector-fold-left
          (lambda (overlaps syn-low syn-high)
            (let ((synapse (synapses-search synapses syn-low syn-high)))
              (if (and synapse (fx>=? (syn-perm synapse) threshold))
                  (add1 overlaps)
                  overlaps)))
            0
          syn-lows syn-highs))
      permanences)))

;; === Column Pooler algorithm ===
                                                                                            ;
(define  compute                         ;; ...
  ;; run one time step of the column pooler algorithm
  (case-lambda
#;> [(cp feedforward-input learn)        ;; CP {CellX} Boolean ->
      (compute cp feedforward-input (list) (list) learn (list)) ]
                                                                                            ;
#;> [(cp feedforward-input               ;; CP {CellX} {{CellX}} {CellX} Boolean ->
        lateral-inputs feedforward-growth-candidates learn)
      (compute cp feedforward-input lateral-inputs feedforward-growth-candidates learn (list)) ]
                                                                                            ;
#;> [(cp feedforward-input               ;; CP {CellX} {{CellX}} {CellX} Boolean {CellX} ->
        lateral-inputs feedforward-growth-candidates learn predicted-input)
      (when (fx<? (length (cp-distal-permanences cp)) (length lateral-inputs))
        (cp-distal-permanences-set! cp
          (make-list (length lateral-inputs) (make-vector (cp-cell-count cp) (fxvector)))))
      (let ((feedforward-growth-candidates
              (if (null? feedforward-growth-candidates) feedforward-input
                  feedforward-growth-candidates)))
        (if (cp-online-learning cp)
          (cond                          ;; online learning step
            [(and (pair? predicted-input)
                  (fx>? (length predicted-input) (cp-predicted-inhibition-threshold cp)))
              (let ((predicted-active-input
                      (intersect1d feedforward-input predicted-input))
                    (predicted-growth-candidates
                      (intersect1d feedforward-growth-candidates predicted-input)))
                (compute-inference-mode cp predicted-active-input lateral-inputs)
                (compute-learning-mode cp predicted-active-input lateral-inputs 
                                          feedforward-growth-candidates)) ]
            [(not (fx<=? (cp-min-sdr-size cp)
                         (length (cp-active-cells cp))
                         (cp-max-sdr-size cp)))
              ;; if no single representation, try to infer one before attempting to learn
              (compute-inference-mode cp feedforward-input lateral-inputs)
              (compute-learning-mode cp feedforward-input lateral-inputs
                                        feedforward-growth-candidates) ]
            [else 
              ;; no predicted input and single SDR, extend that representation
              (compute-learning-mode cp feedforward-input lateral-inputs
                                        feedforward-growth-candidates) ])
          ;; not online-learning
          (if learn
                                         ;; learning step
            (compute-learning-mode cp feedforward-input lateral-inputs 
                                      feedforward-growth-candidates)
                                         ;; inference step
            (compute-inference-mode cp feedforward-input lateral-inputs)))) ] ))
                                                                                           ;
(define (compute-learning-mode cp        ;; CP {CellX} {{CellX}} {CellX} ->
          feedforward-input lateral-inputs feedforward-growth-candidates)
  ;; in learning mode, maintain prior activity or create random sdr for new object
  (define (new-sdr)                      ;; -> SDR
    ;; produce random sdr with exactly sdr-size bits
    (let loop [(n 0) (xs (list))]
      (cond [(fx=? n (cp-sdr-size cp)) (list-sort fx<? xs)]
            [else
              (let ((try (random (cp-cell-count cp))))
                (if (memv try xs)
                  (loop n xs)
                  (loop (add1 n) (cons try xs))))])))
  (let ((prev-active-cells (cp-active-cells cp)))
    (when (cond
        [ (fx<? (length (cp-active-cells cp)) (cp-min-sdr-size cp))
            ;; not enough previously active cells: create a new sdr and learn on it
            (cp-active-cells-set! cp (new-sdr))
            #t ]
        [ (fx>? (length (cp-active-cells cp)) (cp-max-sdr-size cp))
            ;; union of cells active: don't learn
            #f ]
        [ else #t ] )
      ;; do the actual learning
      (when (positive? (length feedforward-input))
        (learn! (cp-proximal-permanences cp)
          (cp-active-cells cp) feedforward-input feedforward-growth-candidates
          (cp-sample-size-proximal cp) (cp-initial-proximal-permanence cp)
          (cp-syn-perm-proximal-inc cp) (cp-syn-perm-proximal-dec cp))
        (do ( (i 0 (add1 i))             ;; external distal learning
              (lateral-input lateral-inputs (cdr lateral-input)))
            ((fx=? i (length lateral-inputs)))
          (learn! (list-ref (cp-distal-permanences cp) i)
            (cp-active-cells cp) (car lateral-input) (car lateral-input)
            (cp-sample-size-distal cp) (cp-initial-distal-permanence cp)
            (cp-syn-perm-distal-inc cp) (cp-syn-perm-distal-dec cp)))
        (learn! (cp-internal-distal-permanences cp)
          (cp-active-cells cp) prev-active-cells prev-active-cells
            (cp-sample-size-distal cp) (cp-initial-distal-permanence cp)
            (cp-syn-perm-distal-inc cp) (cp-syn-perm-distal-dec cp))))))
                                                                                            ;
(define (compute-inference-mode cp       ;; CP {CellX} {CellX} ->
          feedforward-input lateral-inputs)
  ;; in inference mode, if there is some feedforward activity, recognize previously
  ;; known objects and activate a subset of the cells with feedforward support;
  ;; otherwise use lateral activity to activate a subset of the previous active cells
  (define (extend-active active          ;; {CellX} {CellX} {Nat} Boolean -> {CellX}
            cells num-active-segs-for-cells allow-zero)
    ;; extend active with cells in order of activation until sdrSizequota is reached
    (let loop [ (ttop (apply fxmax num-active-segs-for-cells))
                (active active) ]
      (if (and (if allow-zero (not (negative? ttop)) (positive? ttop))
               (fx<? (length active) (cp-sdr-size cp)))
        (loop (fx- ttop 1)
              (append
                (fold-left               ;; filter cells with activation = current ttop
                  (lambda (filtered cellx nasfc)
                    (if (fx=? nasfc ttop)
                        (cons cellx filtered)
                        filtered))
                  (list)
                  cells num-active-segs-for-cells)
                active))
        (list-sort fx<? active))))
  (let* ( (prev-active-cells (cp-active-cells cp))
          (overlaps (right-vec-sum-at-nzgte-thresh (cp-proximal-permanences cp)
                        feedforward-input (cp-connected-permanence-proximal cp)))
          (feedforward-supported-cells
            (indices-where fx>=? overlaps (cp-min-threshold-proximal cp)))
          (num-active-segments-by-cell (make-vector (cp-cell-count cp) 0)))
    (let ((overlaps (right-vec-sum-at-nzgte-thresh (cp-internal-distal-permanences cp)
                      prev-active-cells (cp-connected-permanence-distal cp))))
      (increment-where! num-active-segments-by-cell 
                        fx>=? overlaps (cp-activation-threshold-distal cp)))
    (for-each
      (lambda (lateral-input distal-permanence)
        (let ((overlaps (right-vec-sum-at-nzgte-thresh distal-permanence
                          lateral-input (cp-connected-permanence-distal cp))))
          (increment-where! num-active-segments-by-cell 
                            fx>=? overlaps (cp-activation-threshold-distal cp))))
      lateral-inputs (cp-distal-permanences cp))
    (let ((chosen-cells (list)))
      (unless (zero? (length feedforward-supported-cells))
        (let ((num-active-segs-for-ff-sup-cells
                (list-of (vector-ref num-active-segments-by-cell x)
                         (x in feedforward-supported-cells))))
          (set! chosen-cells (extend-active chosen-cells  ;; exclude cells with 0 active segments
                               feedforward-supported-cells num-active-segs-for-ff-sup-cells #f))))
      (when (and (fx<? (length chosen-cells) (cp-sdr-size cp)) (cp-use-inertia cp))
        (let* ( (prev-cells (setdiff1d prev-active-cells chosen-cells))
                (inertial-cap (int<- (* (length prev-cells) (cp-inertia-factor cp)))))
          (when (positive? inertial-cap)
            (let* ( (num-active-segs-for-prev-cells 
                      (list-of (vector-ref num-active-segments-by-cell x)
                               (x in prev-cells)))
                    (sort-indices (list-sort (lambda (x y)
                          (fx>? (list-ref num-active-segs-for-prev-cells x)
                                (list-ref num-active-segs-for-prev-cells y)))
                        (indexes num-active-segs-for-prev-cells)))
                    (prev-cells (take inertial-cap
                        (list-of (list-ref prev-cells x) (x in sort-indices))))
                    (num-active-segs-for-prev-cells (take inertial-cap
                        (list-of (list-ref num-active-segs-for-prev-cells x) 
                                 (x in sort-indices)))))
              (set! chosen-cells (extend-active chosen-cells
                                   prev-cells num-active-segs-for-prev-cells #t))))))
      (let ((discrepancy (fx- (cp-sdr-size cp) (length chosen-cells))))
        (when (positive? discrepancy)
          (let* ( (rem-ff-cells (setdiff1d feedforward-supported-cells chosen-cells))
                  (n (fxdiv (fx* (length rem-ff-cells) discrepancy) (cp-sdr-size cp)))
                  (n (fxmin (fxmax n discrepancy) (length rem-ff-cells))))
            (if (fx>? (length rem-ff-cells) n)
              (let ((selected (vector->list (vector-sample (list->vector rem-ff-cells) n))))
                (set! chosen-cells (append selected chosen-cells)))
              (set! chosen-cells (append rem-ff-cells chosen-cells))))))
      (cp-active-cells-set! cp (unique eqv? (list-sort fx<? chosen-cells))))))
                                                                                            ;
(define (num-connected-proximal-synapses ;; CP {CellX} -> Nat
          cp cells)
  ;; produce total count of synapses above threshold of cells
  (let ((threshold (cp-connected-permanence-proximal cp)))
    (fold-left
      (lambda (total cellx)
        (vector-fold-left
          (lambda (total synapse)
            (if (fx>=? (syn-perm synapse) threshold)
              (add1 total)
              total))
          total
          (vector-ref (cp-proximal-permanences cp) cellx)))
      0
      cells)))
                                                                                            ;
(define (number-of-proximal-synapses     ;; CP {CellX} -> Nat
          cp cells)
  ;; produce total count of synapses for cells
  (fold-left
    (lambda (total cellx)
      (fx+ total
         (synapses-length (vector-ref (cp-proximal-permanences cp) cellx))))
    0
    cells))
                                                                                            ;
(define (reset cp)                       ;; CP ->
  ;; when learning this signifies we are to learn a unique new object
  (cp-active-cells-set! cp (list)))
                                                                                            ;
(define (learn! permanences              ;; Permanences {CellX} {InputX} {InputX} Nat Perm Perm Perm ->
          active-cells active-input growth-candidate-input
          sample-size initial-permanence permanence-increment permanence-decrement)
  ;; for each active cell, reinforce active synapses, punish inactive synapses,
  ;; and grow new synapses to input bits that the cell isn't already connected to
  (adapt-synapses! permanences active-cells active-input permanence-increment permanence-decrement)
  (if (negative? sample-size)
    (add-synapses! permanences active-cells active-input initial-permanence)
    (let ((max-new-by-cell (max-new-by-cell permanences active-cells active-input sample-size)))
      (grow-synapses! permanences active-cells growth-candidate-input max-new-by-cell initial-permanence))))

;; === Vector and list utilities ===
                                                                                            ;
(define (indices-where pred vec value)   ;; (X Y -> Bool) (vectorof X) Y -> {Nat}
  ;; produce sorted list of indices of vec elements for which (pred element value)
  (do [ (vx (fx- (vector-length vec) 1) (fx- vx 1))
        (result (list) (if (pred (vector-ref vec vx) value)
                           (cons vx result)
                           result)) ]
      ((negative? vx) result)))
                                                                                            ;
(define (increment-where! cs pred xs y)  ;; (Vectorof Nat) (X Y -> Boolean) (Vectorof X) Y ->
  ;; mutates cs, incrementing elements corresponding to xs for which (pred x y)
  (do [ (vx 0 (add1 vx)) ]
      ((fx=? vx (vector-length xs)))
      (when (pred (vector-ref xs vx) y)
        (vector-set! cs vx (add1 (vector-ref cs vx))))))
                                                                                            ;
(define-syntax fold-of                   ;; Op Base Expr Clause ...
  ;; List comprehensions from the Programming Praxis standard prelude
  ;; (see https://programmingpraxis.com/contents/standard-prelude/ for details)
  ;; Only (list-of (var in list)) is used
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

)
