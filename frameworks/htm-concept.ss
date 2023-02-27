;; © 2019 Roger Turner <https://github.com/rogerturner/HTM-scheme/issues/new/choose>
;; SPDX-License-Identifier: AGPL-3.0-or-later  (see Notices below)

#| HTM-scheme Concept (Introductory notes, data structures and core functions) 

Scheme[1] translation of HTM algorithms and experiments used in Numenta research papers[2].

The objective is to run HTM neuroscience experiments with biologically plausible parameters
in minimal memory, reasonable time, on personal hardware, without inscrutable dependencies.

Memory requirement is ~10 bytes/synapse, with simulations of up to 1024 (cortical) columns,
16384 cells/cc (128 minicols of 128 cells in up to 8 (sub-)layers, or equivalent multiple).

Code is R6RS[3] with minor ChezScheme[4] extensions and doesn't use any external libraries.

Algorithms are translated from Numenta "htmresearch" code (using different data structures)
wrapped in framework code to model biological features (eg layer4.ss and coordinates.ss for
sublayers of different cell types and "hexagonal" minicolumn lattice topology).

The key data structures[5] are Segment (dendrite+synapses) and AxonTree (pre-post mapping).
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
using corresponding functions, variable names, and organization (not idiomatic Scheme). [8]
Code is generally "plain Scheme" with no use of continuations or syntax extensions†. Fixnum
operations are used wherever possible, and mutating / proper-list assuming versions of some
standard procedures are included with the utility functions in frameworks/htm-prelude.ss[9]

  † htm-prelude has syntax extensions for examples, "smoke tests", and list comprehensions.

Notes:

  [1] https://en.wikipedia.org/wiki/Scheme_(programming_language)
      "The greatest single programming language ever designed" [Alan Kay]
      "Lisp's parentheses are the bumps on the top of Lego" [Paul Graham]
      "I intend this but for a Scheme of a larger Design" [John Woodward]
      "car and cdr are the only honest function names"  [Citation needed]
      
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
      
  [9] "The optimization will continue until morale improves."
      [Venkatesh Rao]
      
  [A] "Don't worry, you don't have to start your code from scratch."
      https://phdcomics.com/comics/archive.php?comicid=1689 (note "Emergency Button")

Code formatting and idioms:

  Indentation facilitates using an editor "Fold All" view for a file overview - blank lines
  with right margin ";" provide foldable vertical spacing.
  
  Function and parameter names generally follow Numenta code (transformed to "kebab-case").
  Libraries export plain (internal) names: using modules may prefix or rename on importing.
  
  Some function names in Numenta/numpy code are retained (eg. union1d), although the Scheme
  version is specialized for use in HTM-scheme (and could reasonably be renamed SDR-union).

  Function definitions are usually commented only with their type and one-line description,
  and optionally with '(example: (fn args) => result )' form(s), which are checked on load.
  For core algorithms it may be useful to view Scheme and corresponding Numenta Python code
  side-by-side to see Numenta comments for fuller descriptions of functions and parameters.
  
  Comments may include cross-references by [tagref](https://github.com/stepchowfun/tagref):
  "[tag: name] note" can be referenced from "[ref: name] usage" when usage depends on note.

  Some libraries have checked examples and "smoke tests", executed on loading the library.
  (Not characterisation tests: may be useful as examples of usage)
  
  Scheme treats any non-#f value as true, and the `(and ...)` form short circuits, so code
  like `(and connect? (connect? arg ...))` evaluates to #f if connect? is #f, or the value
  produced by applying it if it is a procedure.
  
  "One-armed" `if` is not used, but a simple consequent may follow on the same line as the
  test, with the alternative doubly indented. (Consequent on new line is singly indented.)

  Some key parameters (eg number of minicolumns/cortical column) are set in a "parameters"
  file, which is typically created using command-line arguments by the compilation script.
  
  Record types can be instantiated using the key-word-args procedure to convert parameters
  specified as an unordered list of (key . value) pairs to the constructor arguments, with
  default values for unspecified parameters. Fields which are functions of other arguments
  can be constructed by appending to the arg list and reapplying the constructor (enabling
  them to be immutable), so the protocol clause (function to create a record) may be like:
    (protocol
      (lambda (new)
        (lambda (kwargs)
          (let* ([record (apply new (key-word-args kwargs defaults))]
                 [dependent-parameter (function-of record)] ...
                 [kwargs (append dependent-parameters kwargs)])
            (apply new (key-word-args kwargs defaults))))))
  (Note: applying new twice may be outwith the r6rs standard, but Chez Scheme accepts it.)

Types:
  Pair, Number, Integer, Boolean, Fixnum, (Listof X), (Vector X->Y) ... = Scheme types
  (X Y -> Z)   function with argument types X, Y and result type Z
  {X}          abbreviation for (Listof X)
  Nat          natural number (including zero) (Scheme Fixnum or exact Integer)
  Perm         Nat permanence: 0-255 interpreted as 0.0-1.0 (resolution ~0.005)
  CCX          Nat cortical column index (typically 0-1023)
  LayerX       Nat index identifying layer or cell population (typically 0-7)
  CellX        Nat index of cell in layer (typically 0-2047)
  ColX         Nat minicolumn index of cell: cellx div cells/minicolumn for this layer
  Source       Nat presynaptic cell identifier: ccx <+ layerx <+ cellx
  Synapse      Nat HTM synapse: source <+ perm
  Synapses     Bytevector of Synapse: 32-bit elements, sorted
  Segment      Record with CCX, CellX, Synapses, and overlap counts (40 + 4*nSynapses bytes)
  SegX         Nat 24-bit index of basal/apical segment within layer
  SegXvec      Bytevector of SegX (extendable: bytes 0-1 are index of last used element)
  ColXmask     Bytevector, 1 bit per minicolumn in layer
  Target       Bytevector combining ColXmask and SegXvec
  AxonTree     [(Hashtable Source->Target) . (Vector SegX->Segment)]
  Layer        Record with algorithm parameters, etc
  Macrocolumn  structure of Layers with interconnections (cortical column)
  Patch        multiple Macrocolumns

  |#

  #!chezscheme

(library (frameworks htm-concept)
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
  neighbours-mask
  segx-bytes
  segxv-last
  segxv-ref
  segxv-push
  segxv-map
  segxv-remove!
  make-at
  at-seg-ref
  at-make-seg
  at-free-seg
  at-target
  at-update
  sort-unique!
  sort-unique-by!
  )
                                                                                            ;
(import 
  (chezscheme)
  (parameters)
  (frameworks htm-prelude)
  (frameworks sorting)
  (frameworks coordinates))
                                                                                            ;
  (implicit-exports #f)
  
  ;; Convenience abbreviations
  (alias fxasl  fxarithmetic-shift-left)
  (alias fxasr  fxarithmetic-shift-right)
  (alias bvu32@ bytevector-u32-native-ref)
  (alias bvu32! bytevector-u32-native-set!)
  (define native (native-endianness))
  
  #| Parameters: minicolumns/macrocolumn, deltille-topology?, perm-bits,
                 source-bits, ccx-bits, layer-bits, cellx-bits, and segx-bits
     are imported from parameters.ss)  |#

;; --- Permanences, Sources, Synapses ---
                                                                                            ;
  ;; Synapse is Source<+Perm, where Source is CCX<+LayerX<+CellX
                                                                                            ;
(define min-perm                         ;; Nat
  0)
                                                                                            ;
(define max-perm                         ;; Nat
  (- (expt 2 perm-bits) 1))
  (example: max-perm => #xff )           ;; perm-bits 8
                                                                                            ;
(define (clip-max perm)                  ;; Perm -> Perm
  (fxmin max-perm perm))
                                                                                            ;
(define (clip-min perm)                  ;; Perm -> Perm
  (fxmax min-perm perm))
                                                                                            ;
(define (perm x)                         ;; Number[0.0-1.0] -> Perm
  ;; 0.0 -> 0, but don't round down non-zero x
  (example: (perm 0.0 ) => 0 )
  (example: (perm .001) => 1 )
  (example: (perm .005) => 1 )
  (example: (perm .006) => 2 )
  (example: (perm 0.5 ) => 128 )
  (example: (perm 1.0 ) => 255 )
  (if (zero? x)  0
      (clip-max (fxmax 1 (int<- (* x max-perm))))))
                                                                                            ;
(define source-mask                      ;; Fixnum
  ;; mask for source in synapse
  (fxasl (- (expt 2 source-bits) 1) perm-bits))
  (example: source-mask => #xffffff00 )  ;; source-bits 24
                                                                                            ;
(define (make-source ccx layerx cellx)   ;; CCX LayerX CellX -> Source
  (example: (make-source 1 1 1) => #x004801 )
  (fx+ (fxasl (fx+ (fxasl ccx layer-bits) layerx) cellx-bits) cellx))
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
  (example: (make-syn 7 128 #;(perm 0.5)) => #x00000780 )
  (fx+ (fxasl source perm-bits) perm))
                                                                                            ;
(define (syn-source synapse)             ;; Synapse -> Source
  (fxasr synapse perm-bits))
                                                                                            ;
(define (syn-perm synapse)               ;; Synapse -> Perm
  (fxand synapse max-perm))

;; --- Segments and Synapses ---
                                                                                            ;
(define-record-type seg                  ;; Segment: SegX+CellX, Synapses, Overlaps
  (fields
    (immutable xsource xsource)          ;; Fixnum: segx<+source that this is a segment of
    (mutable synapses)                   ;; Bytevector[u32]: the segment's Synapses
    (mutable overlap))                   ;; Fixnum: overlaps (see calculate-segment-activity)
  (sealed #t) (opaque #t) (nongenerative seg)
                                                                                            ;
(protocol #;(make-seg segx cellx)        ;; SegX CellX -> Segment
  ;; produce a new segment
  (lambda (new)
    (lambda (segx cellx)
      (new (fx+ (fxasl segx source-bits) cellx) (make-synapses 0) 0)))))
                                                                                            ;
(define (seg-cellx seg)                  ;; Segment -> CellX
  (source-cellx (xsource seg)))
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
(define (synapses-set! bv n syn)         ;; Synapses Nat Synapse ->
  (bvu32! bv (fx* n 4) syn))
                                                                                            ;
(define (synapses-length bv)             ;; Synapses -> Nat
  (fxdiv (bytevector-length bv) 4))
                                                                                            ;
(define synapses (lambda source+perm     ;; {Source Perm} -> Synapses
  ;; produce synapses bytevec from alternating source, perm
  (example: (synapses 1 1 2 3) => #vu8(1 1 0 0 3 2 0 0) )
  (let ([synapses (make-synapses (fxdiv (length source+perm) 2))])
    (do ( [s+p source+perm (cddr s+p)]
          [sx 0 (fx+ 4 sx)])
        ((null? s+p) synapses)
      (bvu32! synapses sx (make-syn (car s+p) (cadr s+p)))))))
                                                                                            ;
(define (seg-synapses->list seg)         ;; Segment -> {Synapse}
  (bytevector->uint-list (seg-synapses seg) native 4))
                                                                                            ;
(define (synapses-for-each proc seg)     ;; (Synapse -> ) Segment ->
  ;; apply proc to synapses
  (let* ( [synapses (seg-synapses seg)]
          [n-syns   (synapses-length synapses)])
    (do ([i 0 (fx1+ i)]) ((fx=? i n-syns))
      (proc (synapses-ref synapses i)))))
                                                                                            ;
(define (synapses-map proc seg)          ;; (Synapse -> X) Segment -> {X}
  ;; map synapses by proc
  (map proc
    (bytevector->uint-list (seg-synapses seg) native 4)))
                                                                                            ;
(define (synapses-count pred? segment)   ;; (Synapse -> Boolean) Segment -> Nat
  ;; produce count of synapses for which (pred? synapse) is not #f
  (let* ( [synapses (seg-synapses segment)]
          [n-syns   (synapses-length synapses)])
    (do ( [i 0 (fx1+ i)]
          [n 0 (if (pred? (synapses-ref synapses i)) (fx1+ n) n)])
        ((fx=? i n-syns) n))))
                                                                                            ;
(define (synapses-search source segment) ;; Source Segment -> Synapse | #f
  ;; binary search for synapse [no benefit from manual inlining, unroll?]
  (let* ( [synapses (seg-synapses segment)]
          [target   (fxasl source perm-bits)])
    (let bisect ([left 0] [right (fx- (bytevector-length synapses) 4)])
      (and (fx<=? left right)
        (let* ( [mid      (fxasl (fxasr (fx+ left right) 3) 2)]
                [synapse  (bvu32@ synapses mid)])
          (cond 
            [ (fx<? synapse target)                     (bisect (fx+ mid 4) right)]
            [ (fx<? target (fxand source-mask synapse)) (bisect left (fx- mid 4)) ]
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
  ;; (24 bit Segx values to allow >32 segs/cell with eg 127 minicols, 16 cells per mc/layer.)
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
  ;; colxmask, last, 2 free (=> Bytevector is 32b for 127 mcol)
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
    (let next-colx ([colxs colxs])
      (cond
        [(null? colxs) mask ]
        [else
          (target-colx-bit! mask (car colxs))
          (next-colx (cdr colxs)) ]))))
                                                                                            ;
(define (target-and? target mask)        ;; Target ColXmask -> Boolean
  ;; produce whether mask intersects target
  (let next-32 ([bvx 0])
    (cond
      [(fx>=? bvx colxmask-length)  #f ]
      [(fxzero? (fxand (bvu32@ target bvx) (bvu32@ mask bvx)))
        (next-32 (fx+ 4 bvx)) ]
      [else  #t ] )))
                                                                                            ;
(define all-mask                         ;; ColXmask
  (make-bytevector colxmask-length #xff))
                                                                                            ;
(define (target-zero? target)            ;; Target -> Boolean
  ;; produce whether mask zero
  (not (target-and? target all-mask)))
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
  ;; replace element at index segxx with segx
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
    (do ( [segxx  segxv-base (fx+ segxx segx-bytes)]
          [result (list)     (cons (proc (segxv-ref target segxx)) result)])
        ((fx>? segxx lastx) result))))
                                                                                            ;
(define (segxv-remove! proc target)      ;; (SegX -> Boolean) Target ->
  ;; remove element of segxv for which proc returns #t
  (let ([lastx (segxv-last target)])
    (let loop ([segxx lastx])
      (cond
        [(fx<? segxx segxv-base) ]
        [(proc (segxv-ref target segxx))
          (bytevector-copy! target (fx+ segxx segx-bytes) target segxx (fx- lastx segxx))
          (segxv-last! target (fx- lastx segx-bytes)) ]
        [else  (loop (fx- segxx segx-bytes)) ]))))

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
(define (at-seg-ref at segx)             ;; AxonTree SegX -> Segment
  ;; produce segment segx of at
  (vector-ref (cdr at) segx))
                                                                                            ;
(define (at-make-seg at cellx)           ;; AxonTree CellX -> Segment
  ;; produce a new segment, adding it to the at's segment table
  ;; use a free entry if available, or extend segment table
  (let* ( [seg-table (cdr at)]
          [segx      (vector-ref seg-table 0)]
          [seg       (make-seg segx cellx)])
    #;(assert (fx<? segx (expt 2 segx-bits)))
    (let ([free-next (vector-ref seg-table segx)])
      (vector-set! seg-table segx seg)
      (if (fxpositive? free-next)  (vector-set! seg-table 0 free-next)
          (let* ( [cur-length (vector-length seg-table)]
                  [seg-table  (vector-extend seg-table 1024)]
                  [new-length (vector-length seg-table)])
            (vector-set-fixnum! seg-table (fx1- new-length) 0)
            (do ([sx (fx- new-length 2) (fx1- sx)]) ((fx<? sx cur-length))
              (vector-set-fixnum! seg-table sx (fx1+ sx)))
            (vector-set-fixnum! seg-table 0 cur-length)
            (set-cdr! at seg-table))))
    seg))
                                                                                            ;
(define (at-free-seg at segx)            ;; AxonTree SegX ->
  ;; add segment table entry for segx to free list
  (let ([seg-table (cdr at)])            ;; trap double-free
    (assert (not (fixnum? (vector-ref seg-table segx))))
    (vector-set-fixnum! seg-table segx (vector-ref seg-table 0))
    (vector-set-fixnum! seg-table 0 segx)))
                                                                                            ;
(define (at-target at source)            ;; AxonTree Source -> Target | #f
  ;; produce target for source, or #f if not found
  (eq-hashtable-ref (car at) source #f))
                                                                                            ;
(define (at-new-segment-cells at cellxs) ;; AxonTree {CellX} -> 
  ;; ?
  #f)
                                                                                            ;
(define (at-update at                    ;; AxonTree ColX {Source} Perm Segment (Synapse -> Synapse) ->
           colx ss perm seg adjust)
  ;; merge synapses made from ss (which is sorted) + perm into seg, update existing synapses
  ;; main learning proc: add "growth candidates" + adjust "reinforce candidates" in one pass
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
                        (eq? seg (vector-ref (cdr at) segx)))
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

;; --- Smoke tests ---
                                                                                            ;
  (let ([seg (make-seg 0 0)])
    (seg-synapses-set! seg (synapses 10 0 11 1 12 2 13 3 24 4 25 5 26 6))
    (for-each (lambda (s)
        (expect [(synapses-search s seg)
                 (case s
                  [10 #xA00] [11 #xB01] [12 #xC02] [13 #xD03]
                  [24 #x1804] [25 #x1905] [26 #x1A06]
                  [else #f])]))
      (iota 30)))
                                                                                            ;
  (let ([seg (make-seg 0 0)])
    (seg-synapses-set! seg (synapses 0 0 2 22 4 44))
    [expect ([seg-synapses seg]  #vu8(0 0 0 0 22 2 0 0          44 4 0 0) )]
    [synapses-merge! (list 2 3) 33 seg]
    [expect ([seg-synapses seg]  #vu8(0 0 0 0 22 2 0 0 33 3 0 0 44 4 0 0) )]
    [synapses-merge! (list 3 4) 33 seg]
    [expect ([seg-synapses seg]  #vu8(0 0 0 0 22 2 0 0 33 3 0 0 44 4 0 0) )] )
                                                                                            ;
  (let ([segxv0 (make-target)])          ;; requires segxv-base 18, segx-bytes 3
    [expect ([bytevector-length segxv0]   24 )
            ([segxv-push 1 segxv0]  #f )
            ([segxv-last segxv0]    18 )
            ([segxv-ref segxv0 18]   1 )
            ((not [segxv-push 2 segxv0])  #f )
            ([segxv-last segxv0]    21 )
            ([segxv-ref segxv0 18]   1 )
            ([segxv-ref segxv0 21]   2 )])
                                                                                            ;
  (let ([at (make-at 99)])
    [expect ((cdr at)  '#(1 0) )]
    (let* ( [seg1 (at-make-seg at 1)]
            [ss   (cdr at)])
      [expect ((vector-take 5 ss)  `#(2 ,seg1 3 4 5) )
              ((vector-ref ss (- (vector-length ss) 1))  0 ) ]
      (let ([seg2 (at-make-seg at 2)])
        [expect ((vector-take 5 ss)  `#(3 ,seg1 ,seg2 4 5) ) ]
                 (at-free-seg at 1)
        [expect ((vector-take 5 ss)  `#(1 3 ,seg2 4 5) ) ])))

;; (not effective on import, but will run when any export used):
;; (eval-when (compile eval load visit revisit)
     (check-examples)
;; )

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
