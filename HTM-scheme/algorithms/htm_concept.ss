;;  HTM-scheme Concept  (C) 2019-2021 Roger Turner.  https://discourse.numenta.org/u/rogert
#|  License: GNU Affero Public License version 3  https://www.gnu.org/licenses/agpl-3.0.txt

Scheme[1] translation of HTM algorithms and experiments used in Numenta research papers[2].
The objective is to run HTM neuroscience experiments with biologically plausible parameters
in minimal memory (eg. ~6 bytes/synapse for 100 cortical columns, 2^20 cells 2^30 synapses)

Code is R6RS[3] with minor ChezScheme[4] extensions and doesn't use any external libraries.
Algorithm implementations are as close as possible to Numenta "htmresearch" code wrapped by
framework code to model biological features (eg. layer4.ss and coordinates.ss for hexagonal
minicolumn lattice with cross cortical column connectivity).

The main data structures are Segment (modelling dendrites) and AxonTree (pre-post mapping).
Segments contain sorted vectors of Synapses (presynaptic Source and Permanence in 32 bits),
cell index and "time-stamped" overlap count. Source can be in any cortical column or layer.
AxonTrees are one-to-many mappings (not trees) of Sources to Segments. They are implemented
as a pair of a Hashtable (with key Source, value a vector of SegXs: 16-bit segment numbers)
with a segment table (Vector SegX->Segment). SDRs are lists of column or cell index number. 
The overall TM flow can be sketched as:

->Context input [all active pre-synaptic cells]                \
---> ... [lookup of Source in AxonTree] -> Segments             \
-----> Segments ... [search for synapse for input]               ) depolarize-cells
-------> Synapse with Permanence ... [save overlap in segment]  /
---------> Active/matching Segments for learning, predictions  /
->Feedforward input, learning cells                            \
---> ... [learning: adjust synapses, punish segments, etc]      ) activate-cells
-----> ... [grow new segments and synapses per connect rules]  /

Core algorithms have been translated from github.com/numenta/htmresearch and ../nupic
*using corresponding functions, names, and organization* (not idiomatic Scheme).

Code formatting and idioms:

  Indentation facilitates using a "Fold All" view (eg. Atom) for a file overview (the lines
  with right margin ; provide foldable vertical spacing; #;> is a comment for indentation).
  
  Function and parameter names generally follow Numenta code (transformed to "kebab-case").
  Libraries export plain (internal) names: using modules may prefix or rename on importing.

  Function definitions are usually commented only with their type and one-line description.
  For core algorithms it may be useful to view Scheme and corresponding Numenta Python code
  side-by-side to see Numenta comments for fuller descriptions of functions and parameters.

  Some libraries have "smoke tests" at the end, executed when the library is loaded.
  (Not characterisation tests: may be useful as examples of usage)

  Records can be instantiated using the key-word-args procedure to convert parameters
  specified as an unordered list of (keyword . value) pairs to the arguments to the
  constructor, with default values for unspecified parameters. Fields which depend on
  other parameters can be constructed by prepending to the arg list and re-applying the
  constructor (enabling them to be immutable), so the protocol clause may be like:
    (protocol
      (lambda (new)
        (lambda (kwargs)
          (let* ([record (apply new (key-word-args kwargs defaults))]
                 [dependent-parameters (some-function-of record)]
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
  ColX         Nat minicolumn index of cell: cellx div cells/minicolumn-layer
  Source       Nat presynaptic cell identifier: ccx<< || layerx<< || cellx
  Synapse      Nat HTM synapse: source<<8 || Perm
  Synapses     Bytevector of Synapse: 32-bit elements, sorted
  Segment      Record with CCX, CellX, Synapses, and overlap counts; 20 + 4*nSynapses bytes
  SegX         Nat 16-bit index of segment in layer
  SegXvec      Bytevector of SegX (extendable: offset of last used element at offset 0)
  AxonTree     [(Hashtable Source->SegXvec) . (Vector SegX->Segment)]
  Layer        Record with algorithm parameters, etc
  Macrocolumn  structure of Layers with interconnections (cortical column)
  Patch        multiple Macrocolumns
                                                                                            ;
Notes:
                                                                                            ;
  [1] https://en.wikipedia.org/wiki/Scheme_(programming_language)
      "The greatest single programming language ever designed" [Alan Kay]
      "Lisp's parentheses are the bumps on the top of Lego" [Paul Graham]
      "car and cdr are the only honest function names"  [citation needed]
      
  [2] https://github.com/numenta/htmpapers

  [3] Dybvig 2009 The Scheme Programming Language 4th Edition (https://www.scheme.com/tspl4/)
      "Kent Dybvig's TSPL is to Scheme what K&R is to C" [Daniel P Friedman]

  [4] https://github.com/cisco/ChezScheme (Apache License 2.0)

"A theory should not attempt to explain all the facts because some of the facts are wrong."
 [Francis Crick]
"Just because you've implemented something doesn't mean you understand it."
 [Brian Cantwell Smith]

  |#

  #!chezscheme

(library (HTM-scheme HTM-scheme algorithms htm_concept)
                                                                                            ;
(export
  cellx-bits
  layer-bits
  ccx-bits
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
  seg-synapses<-list
  synapses-map!
  synapses-for-each
  synapses-search
  synapses-add!
  synapses-update!
  segxv-last
  segxv-ref
  segxv-push
  segxv-memv
  segxv-map
  segxv-remove!
  make-at
  at-make-seg
  at-segxv
  at-seg-ref
  at-seg-set!
  at-update!
  sort-unique!
  sort-unique-n!
  expect
  )
                                                                                            ;
(import 
  (chezscheme)
  (HTM-scheme HTM-scheme algorithms htm_prelude))
                                                                                            ;
  (implicit-exports #f)

  ;; Convenience abbreviations
  (alias fxasl  fxarithmetic-shift-left)
  (alias fxasr  fxarithmetic-shift-right)
  (alias bvu16  bytevector-u16-native-ref)
  (alias bvu16! bytevector-u16-native-set!)
  (alias bvu32  bytevector-u32-native-ref)
  (alias bvu32! bytevector-u32-native-set!)
  (define native (native-endianness))

;; --- Permanence and Synapse values ---
                                                                                            ;
(define perm-bits     8)
                                                                                            ;
(define min-perm      0)
                                                                                            ;
(define max-perm      (- (expt 2 perm-bits) 1))
                                                                                            ;
(define source-shift  perm-bits)
                                                                                            ;
(define source-mask   (fxasl #xFFFFFF source-shift))
                                                                                            ;
(define cellx-bits    15)
(define layer-bits     3)
(define ccx-bits       6)
                                                                                            ;
(define layer-shift   cellx-bits)
(define ccx-shift     (+ layer-bits cellx-bits))
                                                                                            ;
(define (make-source ccx layer cellx)    ;; CCX LayerX CellX -> Source
  (fx+ (fxasl ccx ccx-shift) (fxasl layer layer-shift) cellx))
                                                                                            ;
(define (source-ccx s)                   ;; Source -> CCX
  (fxasr s ccx-shift))
                                                                                            ;
(define (source-layer s)                 ;; Source -> LayerX
  (fxand (fxasr s layer-shift) #x7))
                                                                                            ;
(define (source-cellx s)                 ;; Source -> CellX
  (fxand s #x7FFF))
                                                                                            ;
(define (make-syn source perm)           ;; Source Perm -> Synapse
  (fx+ (fxasl source source-shift) perm))
                                                                                            ;
(define (syn-source synapse)             ;; Synapse -> Source
  (fxasr synapse source-shift))
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

;; --- Segment and Synapses---
                                                                                            ;
(define-record-type seg                  ;; Segment: Synapses, Overlaps, CCX+SegX+CellX
  (fields
    (immutable csc csc)                  ;; Fixnum: ccx/segx/cellx that this is a segment of
    (mutable synapses)                   ;; Bytevector: the segment's Synapses
    (mutable overlap)                    ;; Fixnum: overlaps (see calculate-segment-activity)
    )
  (sealed #t) (opaque #t) (nongenerative seg)
(protocol #;(make-seg ccx segx cellx)    ;; CCX SegX CellX -> Segment
  ;; produce a new segment
  (lambda (new)
    (lambda (ccx segx cellx)
      (new (fx+ (fxasl ccx 32) (fxasl segx 16) cellx) (make-synapses 0) 0)))))
                                                                                            ;
(define (seg-cellx seg)                  ;; Segment -> CellX
  (fxand (csc seg) #x7FFF))
                                                                                            ;
(define (seg-segx seg)                   ;; Segment -> SegX
  (fxand (fxasr (csc seg) 16) #xFFFF))
                                                                                            ;
(define (seg-ccx seg)                    ;; Segment -> CCX
  (fxand (fxasr (csc seg) 32) #xFF))
                                                                                            ;
(define (make-synapses n)                ;; Nat -> Synapses
  (make-bytevector (fx* n 4)))
                                                                                            ;
(define (synapses-ref bv n)              ;; Synapses Nat -> Synapse
  (bvu32 bv (fx* n 4)))
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
(define (seg-synapses<-list seg ss)      ;; Segment {Synapse} ->
  ;; used to add or remove synapses (see also synapses-merge!)
  (seg-synapses-set! seg (uint-list->bytevector ss native 4)))
                                                                                            ;
(define (synapses-map! proc syns)        ;; (Synapse -> Synapse) Synapses ->
  ;; update syns by applying proc
  (do ([i (fx1- (synapses-length syns)) (fx1- i)])
      ((fxnegative? i))
    (synapses-set! syns i (proc (synapses-ref syns i)))))
                                                                                            ;
(define (synapses-for-each proc syns)    ;; (Synapse -> ) Synapses ->
  ;; apply proc to syns
  (do ([i (fx1- (synapses-length syns)) (fx1- i)])
      ((fxnegative? i))
    (proc (synapses-ref syns i))))
                                                                                            ;
(define (synapses-search source syns)    ;; Source Synapses -> Synapse | #f
  ;; binary search for synapse [no benefit from manual inlining, unroll?]
  (let ([target (fxasl source source-shift)])
    (let search ([left 0] [right (fx- (bytevector-length syns) 4)])
      (and (fx<=? left right)
        (let* ( [mid      (fxasl (fxasr (fx+ left right) 3) 2)]
                [synapse  (bvu32 syns mid)])
          (cond 
            [ (fx<? synapse target)                     (search (fx+ mid 4) right)]
            [ (fx<? target (fxand source-mask synapse)) (search left (fx- mid 4)) ]
            [ else synapse ]))))))
                                                                                            ;
(define (synapses-merge! ss nss syns)    ;; {Synapse} Nat Synapses -> Synapses
  ;; merge ss (length nss) into syns, omitting duplicates
  (let* ( [nsyns (bytevector-length syns)]
          [nout  (fx+ nsyns (fx* 4 nss))]
          [out   (make-bytevector nout)])
    (let loop ([ss ss] [synx 0] [outx 0] [dup 0])
      (cond
        [(or (null? ss) (fx=? synx nsyns))
          (if (null? ss)
            (bytevector-copy! syns synx out outx (fx- nsyns synx))
            (do ( [outx outx (fx+ outx 4)] ;; copy tail
                  [ss   ss   (cdr ss)])
                ((null? ss))
              (bvu32! out outx (car ss))))
          (if (fxpositive? dup)
            (bytevector-truncate! out (fx- nout dup))
            out) ]
        [else (let ([synapse (bvu32 syns synx)])
            (cond
              [(fx<? synapse (car ss))
                (bvu32! out outx synapse)
                (loop ss (fx+ synx 4) (fx+ outx 4) dup) ]
              [(fx=? (syn-source synapse) (syn-source (car ss)))
                (bvu32! out outx synapse)
                (loop (cdr ss) (fx+ synx 4) (fx+ outx 4) (fx+ dup 4)) ]
              [else 
                (bvu32! out outx (car ss))
                (loop (cdr ss) synx (fx+ outx 4) dup) ])) ]))))
                                                                                            ;
(define (synapses-add! ss nss perm syns) ;; {Source} Nat Perm Synapses -> Synapses
  ;; merge synapses made from ss (which is sorted and of length nss>0) and perm into syns;
  ;; omit duplicates, use binary search of syns to find where to start merging
  (let* ( [l-in  (bytevector-length syns)]
          [l-out (fx+ l-in (fx* 4 nss))]
          [first (car ss)]
          [out   (make-bytevector l-out)])
    (define (merge inx outx dup)
      (let merge ([inx inx] [outx outx] [dup dup] [ss (cdr ss)])
        (cond
          [(and (pair? ss) (fx<? inx l-in)) 
            (let ([synapse (bvu32 syns inx)]
                  [next    (car ss)])
              (cond
                [(fx<? (syn-source synapse) next)
                  (bvu32! out outx synapse)
                  ;; (merge (fx+ inx 4) (fx+ outx 4) dup ss)
                  (do ([inx (fx+ inx 4) (fx+ inx 4)]
                       [outx (fx+ outx 4) (fx+ outx 4)])
                      ((or (fx=? inx l-in)
                           (let ([synapse (bvu32 syns inx)])
                              (if (fx<? (syn-source synapse) next)
                                (begin (bvu32! out outx synapse) #f)
                                #t)))
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
            (if (fxpositive? dup)
              (bytevector-truncate! out (fx- l-out dup))
              out) ] )))
    (let search ([left 0] [right (fx- l-in 4)])
      (if (fx<=? left right)
        (let* ( [mid   (fxasl (fxasr (fx+ left right) 3) 2)]
                [mid@  (syn-source (bvu32 syns mid))])
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
(define (synapses-update! proc syns)     ;; (Synapse -> Synapse | #f) Synapses -> Synapses
  ;; update syns by applying proc, deleting if proc produces #f; all deleted -> empty bytevec
  (let loop ([synx 0] [l-syns (bytevector-length syns)])
    (if (fx<? synx l-syns)
      (let* ( [old (bvu32 syns synx)]
              [new (proc old)])
        (cond
          [new  (unless (fx=? new old) (bvu32! syns synx new))
                (loop (fx+ synx 4) l-syns) ]
          [else
            (bytevector-copy! syns (fx+ synx 4) syns synx (fx- l-syns synx 4))
            (loop synx (fx- l-syns 4)) ]))
      (bytevector-truncate! syns l-syns))))

;; --- Segment index vectors ---
                                                                                            ;
  ;; SegXvec is Bytevector of 16-bit values, with offset of last value at offset 0, and space
  ;; for more value(s); empty SegXvec is #vu8(0 0 0 0); values index a (Vector SegX->Segment)
                                                                                            ;
(define (make-segxv)                     ;; -> SegXvec
  ;; produce a new empty SegXvec
  (make-bytevector 4 0))
                                                                                            ;
(define (segxv-last segxv)               ;; SegXvec -> Nat [even]
  ;; produce offset of last entry
  (bvu16 segxv 0))
                                                                                            ;
(define (segxv-last! segxv lastx)        ;; SegXvec Nat ->
  ;; store offset of last entry
  (bvu16! segxv 0 lastx))
                                                                                            ;
(define (segxv-ref segxv segxx)          ;; SegXvec Nat -> Nat
  ;; produce entry at offset segxx
  (bvu16 segxv segxx))
                                                                                            ;
(define (segxv-extend segxv)             ;; SegXvec -> SegXvec
  ;; produce copy of segxv with more free space
  (let* ( [len  (bytevector-length segxv)]
          [copy (make-bytevector         ;; (extend to cache-line multiple?)
                  (fx+ len (case len [(4) 20] [(24) 96] [(56) 64] [else 256])))])
    (bytevector-copy! segxv 0 copy 0 len)
    copy))
                                                                                            ;
(define (segxv-push segx segxv)          ;; SegX SegXvec -> SegXvec | #f
  ;; push segx onto segxv: produce new segxv if extended, otherwise #f
  (let ([next-segxx (fx+ 2 (segxv-last segxv))])
    (bytevector-u16-native-set! segxv 0 next-segxx)
    (bytevector-u16-native-set! segxv next-segxx segx)
    (if (fx<? (fx+ 2 next-segxx) (bytevector-length segxv))  #f
        (segxv-extend segxv))))
                                                                                            ;
(define (segxv-memv segx segxv)          ;; SegX SegXvec -> Nat | #f
  ;; produce offset of segx in segxv or #f if not found
  (do ([segxx (segxv-last segxv) (fx- segxx 2)])
      ((or (fxzero? segxx) (fx=? segx (segxv-ref segxv segxx)))
        (if (fxzero? segxx)  #f  segxx))))
                                                                                            ;
(define (segxv-map proc segxv)           ;; (SegX -> X) SegXvec -> {X}
  ;; produce list by applying proc to indices in segxv; last in segxv becomes head of list
  (let ([lastx (segxv-last segxv)])
    (do ( [segxx   2     (fx+ segxx 2)]
          [result (list) (cons (proc (segxv-ref segxv segxx)) result)])
        ((fx>? segxx lastx) result))))
                                                                                            ;
(define (segxv-remove! proc segxv)       ;; (SegX -> Boolean) SegXvec ->
  ;; remove element of segxv for which proc returns #t
  (let ([segxx-last (segxv-last segxv)])
    (let loop ([segxx segxx-last])
      (cond
        [ (fxzero? segxx) ]
        [ (proc (segxv-ref segxv segxx))
            (bytevector-copy! segxv (fx+ segxx 2) segxv segxx (fx- segxx-last segxx))
            (segxv-last! segxv (fx- segxx-last 2)) ]
        [ else  (loop (fx- segxx 2)) ]))))
                                                                                            ;

;; --- AxonTree ---
                                                                                            ;
  ;; AxonTree is [(Hashtable Source->SegXvec) . (Vector SegX->Segment)]
  ;; vector has free space, can be extended: element 0 is the number of entries in use
                                                                                            ;
(define (make-at n-source)               ;; Nat -> AxonTree
  ;; produce axon tree with initial allocation for n-source sources
  `[,(make-eqv-hashtable n-source) . ,(make-vector 2 0)])
                                                                                            ;
(define (at-make-seg at ccx cellx)       ;; AxonTree CCX CellX -> Segment
  ;; produce a new segment, adding it to the at's segment table
  (let* ( [seg-table (cdr at)]
          [segx      (fx1+ (vector-ref seg-table 0))]
          [seg       (make-seg ccx segx cellx)])
    (vector-set-fixnum! seg-table 0    segx)
    (vector-set!        seg-table segx seg)
    (unless (fx<? (fx1+ segx) (vector-length seg-table))
      (set-cdr! at (vector-extend seg-table)))
    seg))
                                                                                            ;
(define (at-segxv at source)             ;; AxonTree Source -> SegXvec | #f
  ;; produce segx vector for source, or #f if not found
  (hashtable-ref (car at) source #f))
                                                                                            ;
(define (at-seg-ref at segx)             ;; AxonTree SegX -> Segment
  ;; produce element segx of at's segment table
  (vector-ref (cdr at) segx))
                                                                                            ;
(define (at-seg-set! at segx seg)        ;; AxonTree SegX Segment ->
  ;; update element segx of at's segment table
  (vector-set! (cdr at) segx seg))
                                                                                            ;
(define (at-update! at source seg)       ;; AxonTree Source Segment ->
  ;; add source, seg to the at if necessary
  (let* ( [cell   (hashtable-cell (car at) source #f)]
          [segx   (seg-segx seg)]
          [segxv  (cdr cell)])
    (if segxv
      (unless (segxv-memv segx segxv)
        (let ([segxv (segxv-push segx segxv)])
          (when segxv
            (set-cdr! cell segxv))))
      (set-cdr! cell (segxv-push segx (make-segxv))))))

;; --- Sorting ---
  ;; Small-list sorting from chezscheme library, for Fixnums and omitting duplicates
  ;; (could use adapted list-merge-sort! from srfi-132 if lists larger?)
                                                                                            ;
(define (sort-unique! ls)                ;; {Fixnum} -> {Fixnum}
  ;; produce ascending sort of ls without duplicates, reusing pairs
  (sort-unique-n! ls (length ls)))
                                                                                            ;
(define (sort-unique-n! ls n)            ;; {Fixnum} Nat -> {Fixnum}
  ;; produce ascending sort of ls (of length n) without duplicates, reusing pairs
  (if (fx<=? n 1)  ls
      (dolsort! ls n (list '()))))
                                                                                            ;
(define (dolsort! ls n loc)              ;; {Fixnum} Nat {{}} -> {Fixnum}
  ;; n is length of ls, loc is a temp for dolmerge!
  (if (fx=? n 1)  (begin (set-cdr! ls '()) ls)
      (let ([i (fxsrl n 1)])
        (let ([tail (list-tail ls i)])
          (dolmerge!
            (dolsort! ls i loc)
            (dolsort! tail (fx- n i) loc)
            loc)))))
                                                                                            ;
(define (dolmerge! ls1 ls2 loc)          ;; {Fixnum} {Fixnum} {{}} -> {Fixnum}
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

  (let ([syns (make-synapses 3)])
    (synapses-set! syns 0 (make-syn 0 00))
    (synapses-set! syns 1 (make-syn 2 22))
    (synapses-set! syns 2 (make-syn 4 44))
    [expect (syns  #vu8(0 0 0 0 22 2 0 0 44 4 0 0) )]
    [expect ([synapses-merge! (list (make-syn 3 33) (make-syn 4 40)) 2 syns]
                   #vu8(0 0 0 0 22 2 0 0 33 3 0 0 44 4 0 0) )]
    [expect ([synapses-add! (list 2 3) 2 33 syns]
                   #vu8(0 0 0 0 22 2 0 0 33 3 0 0 44 4 0 0) )]
    [expect ([synapses-add! (list 3 4) 2 33 syns]
                   #vu8(0 0 0 0 22 2 0 0 33 3 0 0 44 4 0 0) )] )

  (let* ( [segxv0 #vu8(0 0 0 0)]
          [segxv  [segxv-extend segxv0]])
    [expect ([eq? segxv segxv0]           #f )
            ([bytevector-length segxv]    24 )]
    [bytevector-truncate! segxv 4]
    [expect ([bytevector=? segxv segxv0]  #t )])
      
  (let* ( [segxv0 #vu8(0 0 0 0)]
          [segxv1 [segxv-push 1 segxv0]])
    [expect ([eq? segxv1 segxv0]          #f )
            ([bytevector-length segxv1]   24 )
            ([bvu16 segxv1 0]  2 )
            ([bvu16 segxv1 2]  1 )]
    [expect ([segxv-push 2 segxv1]  #f )
            ([bytevector-length segxv1]   24 )
            ([bvu16 segxv1 0]  4 )
            ([bvu16 segxv1 2]  1 )
            ([bvu16 segxv1 4]  2 )])
            
  [expect
    ([sort-unique! '()]         '() )
    ([sort-unique! '(1 2 3)]    '(1 2 3) )
    ([sort-unique! '(2 3 1)]    '(1 2 3) )
    ([sort-unique! '(1 2 2 3)]  '(1 2 3) )
    ([sort-unique! '(0 1 3 5 7 9 0 2 4 6 8 9 9 9)]  '(0 1 2 3 4 5 6 7 8 9) )]

)