;;  HTM-scheme   Copyright 2019-2022 Roger Turner.   https://discourse.numenta.org/u/rogert
#|  License: GNU Affero Public License version 3  https://www.gnu.org/licenses/agpl-3.0.txt

Scheme[1] translation of HTM algorithms and experiments used in Numenta research papers[2].
The objective is to run HTM neuroscience experiments with biologically plausible parameters
in minimal memory (eg ~10 bytes/synapse for 100 cortical columns, 2^20 cells 2^30 synapses)

Code is R6RS[3] with minor ChezScheme[4] extensions and doesn't use any external libraries.
Algorithms are translated from Numenta "htmresearch" code (using different data structures)
wrapped in framework code to model biological features (eg layer4.ss and coordinates.ss for
sublayers of different cell types and hexagonal minicolumn lattice topology).

The main data structures[5] are Segment (dendrite+synapses) and AxonTree (pre-post mapping)

These structures are direct implementations of the HTM neuron model[6] and connectivity[7]:
Segments contain sorted vectors of Synapses (presynaptic Source and Permanence in 32 bits),
cell index and "time-stamped" overlap count. Source can be in any cortical column or layer.
AxonTrees are one-to-many mappings (not trees) of Sources to Segments. They are implemented
as a Hashtable with key Source and value a vector of 24-bit segment numbers, with a segment
table mapping segment number to Segment. Basal and apical segments have separate AxonTrees.

SDRs are sorted Lists of column/cell index numbers. The overall TM flow can be sketched as:

Context input [all active cells] -> Sources SDR             \
-> [lookup each Source in AxonTree] -> Segments              \
---> [search each Segment for Synapse] -> Permanence          } depolarize-cells
-----> [save active/matching overlap count in segment]       /
-------> active/matching Segments, predicted cells SDR      /
Feedforward input [active minicolumns SDR: compare with..   \
-> ..active/matching/predicted] -> bursting, learning cells] \
---> [learning: increment/decrement Permanences..             } activate-cells
-----> ..punish Segments, add Synapses, grow new Segments,   /
-------> .. update AxonTree, -> new active cells SDR]       /

Core algorithms (SP, TM, ATTM, CP) have been translated from numenta/htmresearch and /nupic
using corresponding functions, variable names, and organization (**not idiomatic Scheme**).
Code is generally "plain Scheme", with no use of continuations or syntax extensions. Fixnum
operations are used wherever possible, and mutating / proper-list assuming versions of some
standard procedures are included with the utility functions in frameworks/htm-prelude.ss[8]

Code formatting and idioms:

  Indentation facilitates using a "Fold All" view (eg. Atom) for a file overview (the lines
  with right margin ";" provide foldable vertical spacing)
  
  Function and parameter names generally follow Numenta code (transformed to "kebab-case").
  Libraries export plain (internal) names: using modules may prefix or rename on importing.
  
  Some function names in Numenta/numpy code are retained (eg. union1d), although the Scheme
  version is specialized for use in HTM-scheme (and could reasonably be renamed SDR-union).

  Function definitions are usually commented only with their type and one-line description.
  For core algorithms it may be useful to view Scheme and corresponding Numenta Python code
  side-by-side to see Numenta comments for fuller descriptions of functions and parameters.
  
  Comments may include cross-references by [tagref](https://github.com/stepchowfun/tagref):
  "[tag: name] note" can be referenced from "[ref: name] usage" when usage depends on note.

  Some libraries have "smoke tests" at the end, executed every time the library is loaded.
  (Not characterisation tests: may be useful as examples of usage)
  
  Scheme treats any non-#f value as true, and the "(and ...)" form short circuits, so code
  like "(and connect? (connect? arg ...))" evaluates to #f if connect? is #f, or the value
  produced by applying it if it is a procedure.

  Some key parameters (eg number of minicolumns/cortical column) are set in a "parameters"
  file, which is typically created using command-line arguments by the compilation script.
  
  Record types can be instantiated using the key-word-args procedure to convert parameters
  specified as an unordered list of (key . value) pairs to the constructor arguments, with
  default values for unspecified parameters. Fields which are functions of other arguments
  can be constructed by appending to the arg list and reapplying the constructor (enabling
  them to be immutable), so the protocol clause may be like:
    (protocol
      (lambda (new)
        (lambda (kwargs)
          (let* ([record (apply new (key-word-args kwargs defaults))]
                 [dependent-parameter (function-of record)] ...
                 [kwargs (append dependent-parameters kwargs)])
            (apply new (key-word-args kwargs defaults))))))
                                                                                            ;
Types:
  Pair, Number, Integer, Boolean, Fixnum, (Listof X), (Vector X->Y) ... = Scheme types
  (X Y -> Z)   function with argument types X, Y and result type Z
  {X}          abbreviation for (Listof X)
  Nat          natural number (including zero) (Scheme Fixnum or exact Integer)
  Perm         Nat permanence: 0-255 interpreted as 0.0-1.0
  CCX          Nat cortical column index
  LayerX       Nat index identifying layer or cell population
  CellX        Nat index of cell in layer
  ColX         Nat minicolumn index of cell: cellx div cells/minicolumn for this layer
  Source       Nat presynaptic cell identifier: ccx<< || layerx<< || cellx
  Synapse      Nat HTM synapse: source<<8 || Perm
  Synapses     Bytevector of Synapse: 32-bit elements, sorted
  Segment      Record with CCX, CellX, Synapses, and overlap counts (40 + 4*nSynapses bytes)
  SegX         Nat 24-bit index of basal/apical segment within layer
  SegXvec      Bytevector of SegX (extendable: index of last used element at index 0)
  ColXmask     Bytevector, 1 bit per minicolumn in layer
  Target       Bytevector combining ColXmask and SegXvec
  AxonTree     [(Hashtable Source->Target) . (Vector SegX->Segment)]
  Layer        Record with algorithm parameters, etc
  Macrocolumn  structure of Layers with interconnections (cortical column)
  Patch        multiple Macrocolumns
                                                                                            ;
Notes:
                                                                                            ;
  [1] https://en.wikipedia.org/wiki/Scheme_(programming_language)
      "The greatest single programming language ever designed" [Alan Kay]
      "Lisp's parentheses are the bumps on the top of Lego" [Paul Graham]
      "I intend this but for a Scheme of a larger Design" [John Woodward]
      "car and cdr are the only honest function names"  [citation needed]
      
  [2] https://github.com/numenta/htmpapers

  [3] Dybvig 2009 The Scheme Programming Language 4th Edition (https://www.scheme.com/tspl4/)
      "Kent Dybvig's TSPL is to Scheme what K&R is to C" [Daniel P Friedman]

  [4] https://github.com/cisco/ChezScheme (Apache License 2.0)
      (https://github.com/racket/racket/tree/master/racket/src/ChezScheme for Apple Silicon)
  
  [5] "Show me your [code] and conceal your [data structures], and I will be mystified. Show
       me your [data structures], and I won’t usually need your [code], it’ll be obvious."
      [Fred Brooks (paraphrased)]

  [6] Hawkins & Ahmad 2016 Why Neurons Have Thousands of Synapses, A Theory of Sequence Memory
      in Neocortex (https://doi.org/10.3389/fncir.2016.00023)

  [7] Hawkins Ahmad Cui 2017 A Theory of How Columns in the Neocortex Enable Learning the
      Structure of the World (https://doi.org/10.3389/fncir.2017.00081)

  [8] "Just because you've implemented something doesn't mean you understand it."
      [Brian Cantwell Smith]

  |#

  #!chezscheme

(library (HTM-scheme frameworks htm-concept)
                                                                                            ;
(export
  max-perm
  min-perm
  clip-max
  clip-min
  make-syn
  syn-source
  syn-perm
  perm
  make-source
  source-ccx
  source-layer
  source-cellx
  make-seg
  seg-ccx
  seg-segx
  seg-cellx
  seg-overlap
  seg-overlap-set!
  seg-synapses
  seg-synapses-set!
  make-synapses
  synapses-ref
  synapses-set!
  synapses-length
  seg-synapses->list
  synapses-for-each
  synapses-map
  synapses-count
  synapses-search
  synapses-merge!
  synapses-update!
  segxv-base
  target-cols->mask
  target-and?
  target-adjacent?
  segx-bytes
  segxv-last
  segxv-ref
  segxv-push
  segxv-map
  segxv-remove!
  make-at
  at-make-seg
  at-free-seg
  at-target
  at-seg-ref
  at-seg-set!
  at-update
  sort-unique!
  sort-unique-by!
  expect
  )
                                                                                            ;
(import 
  (chezscheme)
  (parameters)
  (HTM-scheme frameworks htm-prelude)
  (HTM-scheme frameworks coordinates))
                                                                                            ;
  (implicit-exports #f)
  
  ;; Convenience abbreviations
  (alias fxasl  fxarithmetic-shift-left)
  (alias fxasr  fxarithmetic-shift-right)
  (alias bvu32@ bytevector-u32-native-ref)
  (alias bvu32! bytevector-u32-native-set!)
  (define native (native-endianness))
  
  #| (minicolumns/macrocolumn deltille-topology?
      perm-bits source-bits ccx-bits layer-bits cellx-bits segx-bits
      are imported from parameters) |#

;; --- Permanence and Synapse values ---
                                                                                            ;
  ;; Synapse is Source<<Perm, where Source is CCX<<LayerX<<CellX
                                                                                            ;
(define min-perm                         ;; Nat
  0)
                                                                                            ;
(define max-perm                         ;; Nat
  (- (expt 2 perm-bits) 1))
                                                                                            ;
(define source-mask                      ;; Fixnum
  ;; mask for source in synapse
  (fxasl (- (expt 2 source-bits) 1) perm-bits))
                                                                                            ;
(define (make-source ccx layer cellx)    ;; CCX LayerX CellX -> Source
  (fx+ (fxasl (fx+ (fxasl ccx layer-bits) layer) cellx-bits) cellx))
                                                                                            ;
(define (source-cellx s)                 ;; Source -> CellX
  (fxbit-field s 0 cellx-bits))
                                                                                            ;
(define (source-layer s)                 ;; Source -> LayerX
  (fxbit-field s cellx-bits (+ cellx-bits layer-bits)))
                                                                                            ;
(define (source-ccx s)                   ;; Source -> CCX
  (fxbit-field s (+ cellx-bits layer-bits) source-bits))
                                                                                            ;
(define (make-syn source perm)           ;; Source Perm -> Synapse
  (fx+ (fxasl source perm-bits) perm))
                                                                                            ;
(define (syn-source synapse)             ;; Synapse -> Source
  (fxasr synapse perm-bits))
                                                                                            ;
(define (syn-perm synapse)               ;; Synapse -> Perm
  (fxand synapse max-perm))
                                                                                            ;
(define (clip-max perm)                  ;; Perm -> Perm
  (fxmin max-perm perm))
                                                                                            ;
(define (clip-min perm)                  ;; Perm -> Perm
  (fxmax min-perm perm))
                                                                                            ;
(define (perm x)                         ;; Number[0.0-1.0] -> Perm
  ;; 0.0 -> 0, but don't round down non-zero x
  (if (zero? x)  0
      (clip-max (fxmax 1 (int<- (* x max-perm))))))

;; --- Segments and Synapses ---
                                                                                            ;
(define-record-type seg                  ;; Segment: Synapses, Overlaps, CCX+SegX+CellX
  (fields
    (immutable xsource xsource)          ;; Fixnum: segx<<source that this is a segment of
    (mutable synapses)                   ;; Bytevector: the segment's Synapses
    (mutable overlap))                   ;; Fixnum: overlaps (see calculate-segment-activity)
  (sealed #t) (opaque #t) (nongenerative seg)
(protocol #;(make-seg ccx segx cellx)    ;; CCX SegX CellX -> Segment
  ;; produce a new segment
  (lambda (new)
    (lambda (ccx segx cellx)
      (new (fx+ (fxasl segx source-bits) (make-source ccx 0 cellx)) (make-synapses 0) 0)))))
                                                                                            ;
(define (seg-cellx seg)                  ;; Segment -> CellX
  (source-cellx (xsource seg)))
                                                                                            ;
(define (seg-ccx seg)                    ;; Segment -> CCX
  (source-ccx (xsource seg)))
                                                                                            ;
(define (seg-segx seg)                   ;; Segment -> SegX
  (fxbit-field (xsource seg) source-bits (+ source-bits segx-bits)))
                                                                                            ;
(define (make-synapses n)                ;; Nat -> Synapses
  (make-bytevector (fx* n 4)))
                                                                                            ;
(define (synapses-ref bv n)              ;; Synapses Nat -> Synapse
  (bvu32@ bv (fx* n 4)))
                                                                                            ;
(define (synapses-set! bv n u32)         ;; Synapses Nat Synapse ->
  (bvu32! bv (fx* n 4) u32))
                                                                                            ;
(define (synapses-length bv)             ;; Synapses -> Nat
  (fxdiv (bytevector-length bv) 4))
                                                                                            ;
(define (seg-synapses->list seg)         ;; Segment -> {Synapse}
  (bytevector->uint-list (seg-synapses seg) native 4))
                                                                                            ;
(define (synapses-for-each proc seg)     ;; (Synapse -> ) Segment ->
  ;; apply proc to synapses
  (let ([synapses (seg-synapses seg)])
    (do ([i (fx1- (synapses-length synapses)) (fx1- i)])
        ((fxnegative? i))
      (proc (synapses-ref synapses i)))))
                                                                                            ;
(define (synapses-map proc seg)          ;; (Synapse -> X) Segment -> {X}
  ;; map synapses by proc
  (map proc
    (bytevector->uint-list (seg-synapses seg) native 4)))
                                                                  ;
(define (synapses-count pred? segment)   ;; (Synapse -> Boolean) Segment -> Nat
  ;; produce count of synapses for which (pred? synapse) is not #f
  (let ([synapses (seg-synapses segment)])
    (do ( [i (fx1- (synapses-length synapses)) (fx1- i)]
          [n 0 (if (pred? (synapses-ref synapses i)) (fx1+ n) n)])
        ((fxnegative? i) n))))
                                                                                            ;
(define (synapses-search source segment) ;; Source Segment -> Synapse | #f
  ;; binary search for synapse [no benefit from manual inlining, unroll?]
  (let* ( [synapses (seg-synapses segment)]
          [target   (fxasl source perm-bits)])
    (let search ([left 0] [right (fx- (bytevector-length synapses) 4)])
      (and (fx<=? left right)
        (let* ( [mid      (fxasl (fxasr (fx+ left right) 3) 2)]
                [synapse  (bvu32@ synapses mid)])
          (cond 
            [ (fx<? synapse target)                     (search (fx+ mid 4) right)]
            [ (fx<? target (fxand source-mask synapse)) (search left (fx- mid 4)) ]
            [ else synapse ]))))))
                                                                                            ;
(define (synapses-merge! ss perm seg)    ;; {Source} Perm Segment ->
  ;; merge synapses made from ss (which is sorted) and perm into seg;
  ;; omit duplicates, use binary search of syns to find where to start merging
  (let* ( [syns  (seg-synapses seg)]
          [l-in  (bytevector-length syns)]
          [l-out (fx+ l-in (fx* 4 (length ss)))]
          [first (car ss)]               ;; [tag: merge-ss-not-null]
          [out   (make-bytevector l-out)])
    (define (merge inx outx dup)
      (let merge ([inx inx] [outx outx] [dup dup] [ss (cdr ss)])
        (cond
          [(and (pair? ss) (fx<? inx l-in)) 
            (let ([synapse (bvu32@ syns inx)]
                  [next    (car ss)])
              (cond
                [(fx<? (syn-source synapse) next)
                  (bvu32! out outx synapse)
                  ;; (merge (fx+ inx 4) (fx+ outx 4) dup ss)
                  (let run ([inx (fx+ inx 4)] [outx (fx+ outx 4)])
                    (if (fx<? inx l-in)
                      (let ([synapse (bvu32@ syns inx)])
                        (cond
                          [(fx<? (syn-source synapse) next)
                            (bvu32! out outx synapse)
                            (run (fx+ inx 4) (fx+ outx 4)) ]
                          [else (merge inx outx dup ss) ]))
                      (merge inx outx dup ss))) ]
                [(fx>? (syn-source synapse) next)
                  (bvu32! out outx (make-syn next perm))
                  (merge inx (fx+ outx 4) dup (cdr ss)) ]
                [else 
                  (bvu32! out outx synapse)
                  (merge (fx+ inx 4) (fx+ outx 4) (fx+ dup 4) (cdr ss)) ])) ]
          [else
            (if (null? ss)
              (bytevector-copy! syns inx out outx (fx- l-in inx))
              (do ( [outx outx (fx+ outx 4)]  ;; add tail of ss
                    [ss   ss   (cdr ss)])
                  ((null? ss))
                (bvu32! out outx (make-syn (car ss) perm))))
            (seg-synapses-set! seg
              (if (fxpositive? dup)  (bytevector-truncate! out (fx- l-out dup))
                  out)) ] )))
    ;; (synapses-merge!):
    (let search ([left 0] [right (fx- l-in 4)])
      (if (fx<=? left right)
        (let* ( [mid   (fxasl (fxasr (fx+ left right) 3) 2)]
                [mid@  (syn-source (bvu32@ syns mid))])
          (cond
            [ (fx<? mid@  first)  (search (fx+ mid 4) right)]
            [ (fx<? first mid@ )  (search left (fx- mid 4)) ]
            [ else (let ([inx (fx+ mid 4)])
                     (bytevector-copy! syns 0 out 0 inx)
                     (merge inx inx 4)) ]))
        (begin                           ;; 0..left < first: copy then merge
          (bytevector-copy! syns 0 out 0 left)
          (bvu32! out left (make-syn first perm))
          (merge left (fx+ left 4) 0))))))
                                                                                            ;
(define (synapses-update! proc segment)  ;; (Synapse -> Synapse | #f) Segment ->
  ;; update synapses of segment by applying proc, deleting synapse if proc produces #f
  (let* ( [syns   (seg-synapses segment)]
          [l-syns (bytevector-length syns)])
    (let loop ([synx 0] [endx l-syns])
      (if (fx<? synx endx)
        (let* ( [old (bvu32@ syns synx)]
                [new (proc old)])
          (cond
            [new  (unless (fx=? new old) (bvu32! syns synx new))
                  (loop (fx+ synx 4) endx) ]
            [else                        ;; deletion: move following down
              (bytevector-copy! syns (fx+ synx 4) syns synx (fx- endx synx 4))
              (loop synx (fx- endx 4)) ]))
        (when (fx<? endx l-syns)
          (seg-synapses-set! segment (bytevector-truncate! syns endx)))))))

;; --- Targets: projected-to columns and segment index vectors ---
                                                                                            ;
  ;; Each Source has a Targets index (extendable Bytevector) with ColXmask and SegXvec parts.
  ;; ColXmask bits are 1 for minicolumns projected to by that source; SegXvec elements are 24
  ;; bit values which index a (Vector SegX->Segment); last SegX is 16 bits before segxv-base.
                                                                                            ;
(define colxmask-length                  ;; Nat
  ;; padded to mod 4 bytes for target-and?
  (* 4 (quotient (+ 3 (quotient (+ minicolumns/macrocolumn 7) 8)) 4)))
                                                                                            ;
(define segxv-base                       ;; Nat
  ;; index in bytevector of start of SegXvec part (after last index)
  (fx+ colxmask-length 2))
                                                                                            ;  
(define segx-bytes (fxdiv segx-bits 8))  ;; Nat
                                                                                            ;  
  (alias segx@ bytevector-u24-ref)
  (alias segx! bytevector-u24-set!)
                                                                                            ;  
(define initial-target-length            ;; Nat
  ;; colxmask, last, 2 free (32b for 127 mcol )
  (fx+ segxv-base (fx* 2 segx-bytes)))
                                                                                            ;
(define (make-target)                    ;; -> Target
  ;; produce a new Target (ColXmask + lastx + 1+ slots)
  (let ([target (make-bytevector initial-target-length 0)])
    (segxv-last! target (fx- segxv-base segx-bytes))
    target))
                                                                                            ;
(define (target-colx-bit! mask colx)     ;; ColXmask ColX ->
  ;; set bit in mask indicating that colx is projected to
  (bytevector-copy-bit! mask colx 1))
                                                                                            ;
(define (target-cols->mask colxs)        ;; {ColX} -> ColXmask
  ;; produce mask from colxs
  (let ([mask (make-bytevector colxmask-length 0)])
    (for-each (lambda (colx)
        (target-colx-bit! mask colx))
      colxs)
    mask))
                                                                                            ;
(define (target-and? target mask)        ;; Target ColXmask -> Boolean
  ;; produce whether mask intersects target
  (let next-32 ([bvx (fx- colxmask-length 4)])
    (cond
      [(fxnegative? bvx)  #f ]
      [(fxzero? (fxand (bvu32@ target bvx) (bvu32@ mask bvx)))
        (next-32 (fx- bvx 4)) ]
      [else  #t ] )))
                                                                                            ;
(define (target-adjacent? target colx)   ;; Target ColX -> Boolean
  ;; whether target's ColXmask overlaps neighbours-mask of colx ("adjacent" includes colx)
  (target-and? target (vector-ref neighbours-mask colx)))
                                                                                            ;
(define neighbours-mask                  ;; (Vector ColX->ColXmask)
  ;; each mask has 1 bits for colx and 6 adjacent cols
  (build-vector minicolumns/macrocolumn (lambda (colx)
      (do ( [n 0 (+ n 1)]
            [mask (make-bytevector colxmask-length 0)
              (if (>= 1 (within-cc-distance2 colx n))
                (bytevector-copy-bit mask n 1)
                mask) ])
          ((= n minicolumns/macrocolumn) mask )))))
                                                                                            ;
(define (segxv-last target)              ;; Target -> Nat
  ;; produce index of last segx element
  (bytevector-u16-ref target (fx- segxv-base 2) native))
                                                                                            ;
(define (segxv-last! target lastx)       ;; Target Nat ->
  ;; store index of last segx element
  (bytevector-u16-set! target (fx- segxv-base 2) lastx native))
                                                                                            ;
(define (segxv-ref target segxx)         ;; Target Nat -> SegX
  ;; produce element at index segxx
  (segx@ target segxx native))
                                                                                            ;
(define (segxv-set! target segxx segx)   ;; Target Nat SegX ->
  ;; replace element at index segxx with u24
  (segx! target segxx segx native))
                                                                                            ;
(define (segxv-extend target)            ;; Target -> Target
  ;; produce copy of target with more free space
  (let* ( [len  (bytevector-length target)]
          [copy (make-bytevector (fx+ len (fx* 5 segx-bytes)))])
    (bytevector-copy! target 0 copy 0 len)
    copy))
                                                                                            ;
(define (segxv-push segx target)         ;; SegX Target -> Target | #f
  ;; push segx onto target produce new target if extended, otherwise #f
  (let ([next-segxx (fx+ segx-bytes (segxv-last target))])
    (segxv-last! target next-segxx)
    (segxv-set!  target next-segxx segx)
    (if (fx<? (fx+ segx-bytes next-segxx) (bytevector-length target))  #f
        (segxv-extend target))))
                                                                                            ;
(define (segxv-memv segx target)         ;; SegX Target -> Nat | #f
  ;; produce index of segx in target or #f if not found
  (let ([lastx (segxv-last target)])
    (let next ([segxx segxv-base])
      (cond
        [(fx>? segxx lastx) #f ]
        [(fx=? segx (segxv-ref target segxx)) segxx ]
        [else (next (fx+ segxx segx-bytes))]))))
                                                                                            ;
(define (segxv-map proc target)          ;; (SegX -> X) Target -> {X}
  ;; produce list by applying proc to segx values in target
  (let ([lastx (segxv-last target)])
    (do ( [segxx segxv-base (fx+ segxx segx-bytes)]
          [result (list)             (cons (proc (segxv-ref target segxx)) result)])
        ((fx>? segxx lastx) result))))
                                                                                            ;
(define (segxv-remove! proc target)      ;; (SegX -> Boolean) Target ->
  ;; remove element of segxv for which proc returns #t
  (let ([lastx (segxv-last target)])
    (let loop ([segxx lastx])
      (cond
        [ (fx<? segxx segxv-base) ]
        [ (proc (segxv-ref target segxx))
            (bytevector-copy! target (fx+ segxx segx-bytes) target segxx (fx- lastx segxx))
            (segxv-last! target (fx- lastx segx-bytes)) ]
        [ else  (loop (fx- segxx segx-bytes)) ]))))

;; --- AxonTrees ---
                                                                                            ;
  ;; AxonTree is [ (Hashtable Source->Targets) . (Vector SegX->(Segment|Nat)) ]
  ;; element 0 of vector is the head of a free list (chain of indexes of free entries)
  ;; eq-hashtable is used (should be eqv- but keys are always Fixnums)
                                                                                            ;
(define (make-at n-source)               ;; Nat -> AxonTree
  ;; produce axon tree with initial allocation for n-source sources
  (cons* (make-eq-hashtable n-source) (vector 1 0)))
                                                                                            ;
(define (at-make-seg at ccx cellx)       ;; AxonTree CCX CellX -> Segment
  ;; produce a new segment, adding it to the at's segment table
  ;; use a free entry if available, or extend segment table
  (let* ( [seg-table (cdr at)]
          [segx      (vector-ref seg-table 0)]
          [seg       (make-seg ccx segx cellx)])
    ;; (assert (fx<? segx (expt 2 segx-bits)))
    (let ([free-next (vector-ref seg-table segx)])
      (vector-set! seg-table segx seg)
      (if (fxpositive? free-next)  (vector-set! seg-table 0 free-next)
          (let* ( [cur-length (vector-length seg-table)]
                  [seg-table  (vector-extend seg-table 1024)]
                  [new-length (vector-length seg-table)])
            (vector-set-fixnum! seg-table (fx1- new-length) 0)
            (do ([sx (fx- new-length segx-bytes) (fx1- sx)])
                ((fx<? sx cur-length) (vector-set-fixnum! seg-table 0 cur-length))
              (vector-set-fixnum! seg-table sx (fx1+ sx)))
            (set-cdr! at seg-table))))
    seg))
                                                                                            ;
(define (at-free-seg at segx)            ;; AxonTree SegX ->
  ;; add segment table entry for segx to free list
  (let ([seg-table (cdr at)])
    (vector-set-fixnum! seg-table segx (vector-ref seg-table 0))
    (vector-set-fixnum! seg-table 0 segx)))
                                                                                            ;
(define (at-target at source)            ;; AxonTree Source -> Target | #f
  ;; produce target for source, or #f if not found
  (eq-hashtable-ref (car at) source #f))
                                                                                            ;
(define (at-seg-ref at segx)             ;; AxonTree SegX -> Segment
  ;; produce element segx of at's segment table
  (vector-ref (cdr at) segx))
                                                                                            ;
(define (at-seg-set! at segx seg)        ;; AxonTree SegX Segment ->
  ;; update element segx of at's segment table
  (vector-set! (cdr at) segx seg))
                                                                                            ;
(define (at-update at                    ;; AxonTree ColX {Source} Perm Segment (Synapse -> Synapse) ->
           colx ss perm seg adjust)
  ;; merge synapses made from ss (which is sorted) and perm into seg, updating existing synapses
  (define (update-at source)             ;; Source ->
    ;; update ColXmask and SegXvec parts of source's Target
    (let* (                              ;; eq-hashtable-cell will create source entry if necessary
        [ht-cell (eq-hashtable-cell (car at) source #f)]
        [target
          (cond [(cdr ht-cell) ]
                [else (let ([target (make-target)])
                        (set-cdr! ht-cell target)
                        target) ])])
      (target-colx-bit! target colx)
      (unless (segxv-memv (seg-segx seg) target)
        (let ([new-target (segxv-push (seg-segx seg) target)])
          (when new-target
            (set-cdr! ht-cell new-target))))))
  ;; (at-update)
  (let* ( [syns  (seg-synapses seg)]
          [l-in  (bytevector-length syns)]
          [l-out (fx+ l-in (fx* 4 (length ss)))]
          [out   (if (null? ss)  syns            ;; update synapses in place?
                     (make-bytevector l-out))])
    (let merge ([inx 0] [outx 0] [ss ss])
      (define source-fence (fx1+ (syn-source source-mask)))
      (let ([synapse (if (fx=? inx l-in)  (make-syn source-fence 0)  (bvu32@ syns inx))]
            [next-s  (if (null? ss)       source-fence  (car ss))])
        (cond
          [(fx<? (syn-source synapse) next-s)    ;; existing synapse: adjust perm
            (let ([new (adjust synapse)])
              (cond
                [new  (bvu32! out outx new)
                      (merge (fx+ inx 4) (fx+ outx 4) ss) ]
                [else                            ;; #f: remove source's reference to segment
                  (let ([target (at-target at (syn-source synapse))])
                    (segxv-remove! (lambda (segx)
                        (eq? seg (at-seg-ref at segx)))
                      target))
                  (merge (fx+ inx 4) outx ss) ])) ]
          [(fx>? (syn-source synapse) next-s)    ;; new source
            (bvu32! out outx (make-syn next-s perm))
            (update-at next-s)
            (merge inx (fx+ outx 4) (cdr ss)) ]
          [(pair? ss)                            ;; duplicate
            (bvu32! out outx (make-syn next-s (fxmax (syn-perm synapse) perm)))
            (merge (fx+ inx 4) (fx+ outx 4) (cdr ss)) ]
          [else                                  ;; finished
            (unless (and (eq? out syns) (fx=? outx l-out))
              (seg-synapses-set! seg
                (if (fx=? outx l-out)  out
                    (bytevector-truncate! out outx)))) ])))))

;; --- Sorting ---
  ;; Small-list sorting from chezscheme library, for Fixnums and omitting duplicates;
  ;; dolsort! unrolled; could adapt list-merge-sort! from srfi-132 if lists larger?
                                                                                            ;
(define sort-unique!                     ;; {Fixnum} [Nat] -> {Fixnum}
  ;; produce ascending sort of ls [of length n] without duplicates, reusing pairs
  (case-lambda
    [(ls)
      (sort-unique! ls (length ls)) ]
    [(ls n)
      ;; (let () (define (sort-thunk)
      (if (fx<=? n 1)  ls
          (dolsort! ls n (cons '() '())))
      ;; ) (vr0! (fx+ (vr0@) n)) (with-thunk-time sort-thunk))
           ]))
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
(define (sort-unique-by! elt< xs)        ;; (X X -> T?) {X} -> {X}
  ;; produce sort of xs by elt<, without duplicates, reusing pairs
  (let ([n (length xs)])
    ;; (define (sort-thunk)
    (if (fx<=? n 1)  xs
        (dosort-by! elt< xs n (cons '() '())))
    ;; ) (vr0! (fx+ (vr0@) n)) (with-thunk-time sort-thunk)
        ))
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

;; --- Smoke tests ---
                                                                                            ;
(define-syntax expect                    ;; ((X ... -> Y) X ...) Y -> Error |
  ;; check that function application(s) to arguments match expected values
  (lambda (x)                            
    (syntax-case x ()                    ;; [expect ([fn args] expected ) ... ]
      [ (_ (expr expected) ...)          ;; expr matches [fn args]
        #'(begin 
            (let ([result expr])         ;; eval expr once, no output if check passes
              #;(when (equal? result expected) (display "."))
              (unless (equal? result expected)
                (for-each display 
                  `("**" expr #\newline 
                    "  expected: " ,expected #\newline 
                    "  returned: " ,result  #\newline))
                (exit))) ...)])))

  [expect ([perm 0.0  ]    0 )
          ([perm 0.001]    1 )
          ([perm 0.005]    1 )
          ([perm 0.006]    2 )
          ([perm 1.0  ]  255 )]

  (let* ( [seg  (make-seg 0 0 0)]
          [syns (make-synapses 7)])
    (seg-synapses-set! seg syns)
    (synapses-set! syns 0 (make-syn 10 0))
    (synapses-set! syns 1 (make-syn 11 1))
    (synapses-set! syns 2 (make-syn 12 2))
    (synapses-set! syns 3 (make-syn 13 3))
    (synapses-set! syns 4 (make-syn 24 4))
    (synapses-set! syns 5 (make-syn 25 5))
    (synapses-set! syns 6 (make-syn 26 6))
    [expect
      ([synapses-search  0 seg] #f)
      ([synapses-search  9 seg] #f)
      ([synapses-search 10 seg] #xA00)
      ([synapses-search 11 seg] #xB01)
      ([synapses-search 12 seg] #xC02)
      ([synapses-search 13 seg] #xD03)
      ([synapses-search 14 seg] #f)
      ([synapses-search 17 seg] #f)
      ([synapses-search 23 seg] #f)
      ([synapses-search 24 seg] #x1804)
      ([synapses-search 25 seg] #x1905)
      ([synapses-search 26 seg] #x1A06)
      ([synapses-search 27 seg] #f)
      ])
    
  (let* ( [seg  (make-seg 0 0 0)]
          [syns (make-synapses 3)])
    (seg-synapses-set! seg syns)
    (synapses-set! syns 0 (make-syn 0 00))
    (synapses-set! syns 1 (make-syn 2 22))
    (synapses-set! syns 2 (make-syn 4 44))
    [expect ([seg-synapses seg]  #vu8(0 0 0 0 22 2 0 0          44 4 0 0) )]
    [synapses-merge! (list 2 3) 33 seg]
    [expect ([seg-synapses seg]  #vu8(0 0 0 0 22 2 0 0 33 3 0 0 44 4 0 0) )]
    [synapses-merge! (list 3 4) 33 seg]
    [expect ([seg-synapses seg]  #vu8(0 0 0 0 22 2 0 0 33 3 0 0 44 4 0 0) )] )

  #;(let ([segxv0 (make-target)])
    [expect ([bytevector-length segxv0]   24 )
            ([segxv-push 1 segxv0]  #f )
            ([bvu16@ segxv0 16]      18 )
            ([bvu16@ segxv0 18]       1 )
            ([segxv-push 2 segxv0]  #f )
            ([bvu16@ segxv0 16]      20 )
            ([bvu16@ segxv0 18]       1 )
            ([bvu16@ segxv0 20]       2 )])
            
  (let ([at (make-at 99)])
    [expect ([cdr at]  '#(1 0) )]
    (let* ( [seg1 [at-make-seg at 1 1]]
            [ss   [cdr at]])
      [expect ([vector-take 5 ss]  `#(2 ,seg1 3 4 5) )
              ([vector-ref ss [- [vector-length ss] 1]]  0 )]
      (let ([seg2 [at-make-seg at 2 2]])
        [expect ([vector-take 5 ss]  `#(3 ,seg1 ,seg2 4 5) )]
        [at-free-seg at 1]
        [expect ([vector-take 5 ss]  `#(1 3 ,seg2 4 5) )])))
            
  [expect
    ([sort-unique! '()]         '() )
    ([sort-unique! '(1 2 3)]    '(1 2 3) )
    ([sort-unique! '(2 3 1) 3]  '(1 2 3) )
    ([sort-unique! '(1 2 2 3)]  '(1 2 3) )
    ([sort-unique! '(0 9 1 1 1 1 7 5 3 0 2 4 6 1 8 9 9 9 9)]  '(0 1 2 3 4 5 6 7 8 9) )]

  (let ()
    (define (s< s1 s2)
      (string<? (symbol->string s1) (symbol->string s2)))
    [expect
      ([sort-unique-by! s< '()]         '() )
      ([sort-unique-by! s< '(a b c)]    '(a b c) )
      ([sort-unique-by! s< '(b c a)]    '(a b c) )
      ([sort-unique-by! s< '(a b b c)]  '(a b c) )
      ([sort-unique-by! s< '(a c e g c c b a d f h h)]  '(a b c d e f g h) )] )

)