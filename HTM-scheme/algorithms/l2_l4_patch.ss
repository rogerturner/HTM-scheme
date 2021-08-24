;;  === HTM-scheme L2L4 Patch algorithm  (C) 2021 Roger Turner. ===
#|  License: GNU Affero Public License version 3 http://www.gnu.org/licenses

    This program is free software: you can redistribute it and/or modify  
    it under the terms of the GNU Affero Public License version 3 as      
    published by the Free Software Foundation.                            
                                                                         
    This program is distributed in the hope that it will be useful,       
    but WITHOUT ANY WARRANTY; without even the implied warranty of        
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                  
    See the GNU Affero Public License for more details.                   

This L2L4 patch is part of the untangling_sequences project.
See htm_concept.ss for type and data structure description and code conventions.
  
  L2L4 cortical column:

  L2 is a "column pooler" layer, L4 is a "apical tiebreak temporal memory" layer
  (with ss4L4, ss4L23, and p4 cell populations)

                            +---------------<- other L2 cc's active cells
                            v                    (from previous time step)
          +-----------------li--+             
          |          L2        id<--+--------> other L2 cc's lateral inputs  
          +--ff---gc--------ac--+   |            (for next time step)
             ^    ^         v       |              
             |    |         +-------+          ac = active cells
             ^    ^         v                  ai = apical input (to p4 only)
          +--ac---lc--------ai--+              bi = basal input
    +---->bi         L4         |              ff = feedforward input
    |     +-------ff------------+              gc = growth candidates  
    |             ^                            id = internal distal input
    |             |                            li = lateral input (per cc)
    ^             ^                            lc = learning cells
  location      feature

  |#
  
  #!chezscheme

(library (HTM-scheme HTM-scheme algorithms l2_l4_patch)
                                                                                            ;
(export
make-patch
  patch-l2s
  patch-l4s
compute
reset
  reset-l4
  reset-l2
  reset-seq
  print-statistics
  (rename
    (l4:get-active-cells           get-active-cells)
    (l4:get-predicted-cells        get-predicted-cells)
    (l4:get-predicted-active-cells get-predicted-active-cells)
    ))
                                                                                            ;
(import
  (except (chezscheme) reset)
  (parameters)
          (HTM-scheme HTM-scheme algorithms htm_prelude)
          (HTM-scheme HTM-scheme algorithms htm_concept)
          (HTM-scheme HTM-scheme math       coordinates)
  (prefix (HTM-scheme HTM-scheme algorithms layer4) l4:)
  (prefix (HTM-scheme HTM-scheme algorithms layer2) l2:))
                                                                                            ;
  (implicit-exports #f)

  ;; source layer coding:
  ;;    0  "location" (p6l4?)
  ;;    1  apical-input (p23)
  ;;    2  ss4l4
  ;;    3  ss4l23
  ;;    4  p4

;; === L2 L4 algorithm ===
                                                                                            ;
(define-record-type patch (fields        ;; Patch
  l2s                                    ;; CColVecOf L2
  l4s                                    ;; CColVecOf L4
  nearby-ccs)                            ;; VectorOf {Nat} [ccs that are near each cc]
(protocol #;(make-patch ...)
  (lambda (new)
    (lambda (ncc nmc l2-c/mc ss4l4-c/mc ss4l23-c/mc p4-c/mc l6-c/mc cp-overrides attm-overrides)
                                                                                            ;
#;> (define (l2-connect?                 ;; Source ColX CCX -> Boolean
              source                     ;; pre-synaptic cell
              post-mcx                   ;; minicolumn within cc of post-synaptic cell
              post-ccx)                  ;; cortical column of post-synaptic cell
      ;; Could this input intersect with dendrites of a cell in this minicolumn?
      (if no-connectivity?
        (let ([pre-kind (source-layer source)])
          (or (fx=? 1 pre-kind) (fx=? 4 pre-kind)))
        (let* (
            [pre-ccx    (source-ccx source)]
            [pre-kind   (source-layer source)]
            [pre-cellx  (source-cellx source)]
            [pre-mcx    (fxdiv pre-cellx (case pre-kind
                                          [(0) l6-c/mc]
                                          [(1) l2-c/mc]
                                          [(2) ss4l4-c/mc]
                                          [(3) ss4l23-c/mc]
                                          [(4) p4-c/mc]))]
            [same-cc    (fx=? pre-ccx post-ccx)]
            [distance2  (if same-cc
                          (within-cc-distance2 pre-mcx post-mcx)
                          (mc-distance2 pre-ccx pre-mcx post-ccx post-mcx))]
            [~same-mc   (fxpositive? distance2)])
          (case pre-kind
            [(0) (and same-cc (fx<=? distance2 9))]     ;; p6l4 axon span at l2
            [(1) (fx<=? distance2 500)]                 ;; p23 axon span at l2
            [(2) (and same-cc (fx<=? distance2 36))]    ;; ss4l4 axon span at l2
            [(3) (and same-cc (fx<=? distance2 64))]    ;; ss4l23 axon span at l2
            [(4) (fx<=? distance2 500)]))))             ;; p4 axon span at l2
                                                                                            ;
#;> (define (l4-connect?                 ;; Source ColX CCX Layer -> Boolean
              source                     ;; pre-synaptic cell
              post-mcx                   ;; minicolumn within cc of post-synaptic cell
              post-ccx                   ;; cortical column of post-synaptic cell
              post-kind)                 ;; cell/segment type
      ;; Could this input intersect with dendrites of a cell in this minicolumn?
      (let ([pre-kind (source-layer source)])
        (if no-connectivity?
          (case post-kind                ;; no 2d layout: segregated pre->post connections
            [(ss4l4)  (fx=? pre-kind 2)] ;; ss4l4->ss4l4
            [(ss4l23) (fx=? pre-kind 0)] ;; location->ss4l23
            [(p4b)    (fx=? pre-kind 0)] ;; location->p4b
            [(p4a)    (fx=? pre-kind 1)]);; p23->p4a
                                         ;; 2d minicols with distance; apical input is p4 only
          (if (and (fx=? pre-kind 1) (not (eq? post-kind 'p4a)))  #f
            (let* (
                [pre-ccx  (source-ccx source)]
                [cellx    (source-cellx source)]
                [pre-mcx  (fxdiv cellx (case pre-kind
                                        [(0) l6-c/mc]
                                        [(1) l2-c/mc]
                                        [(2) ss4l4-c/mc]
                                        [(3) ss4l23-c/mc]
                                        [(4) p4-c/mc]))]
                [same-cc    (fx=? pre-ccx post-ccx)]
                [distance2  (if (fxnegative? post-mcx)  0
                                (if same-cc
                                  (within-cc-distance2 pre-mcx post-mcx)
                                  (mc-distance2 pre-ccx pre-mcx post-ccx post-mcx)))]
                                         ;; exclude inter-minicol connection?
                [~same-mc   (fxpositive? distance2)]
                ;; reach? = axon/dendrite arbor overlap [Izhikevich & Edelman 2008]:
                ;; (1000 1120 1120 500 150) /50 (minicol spacing) ^2 -> (400 500 500 100 9)
                ;; [Schubert etal 2003]: ss4 axons/dendrites are intra-cc only?
                [reach?
                  (case post-kind
                    [(ss4l4)  (case pre-kind
                        [(0) (and same-cc (fx<=? distance2 400))]  ;; p6l4 axons, ss4l4 dendrites in cc
                        [(1) (fx<=? distance2 9)]                  ;; p23 axon span at l4
                        [(2) (and #;same-cc (fx<=? distance2 500))]  ;; ss4l4 axons within cc
                        [(3) (and same-cc (fx<=? distance2 100))]  ;; ss4l23 axon span at l4
                        [(4) (fx<=? distance2 9)]) ]               ;; p4 axon span
                    [(ss4l23 p4b)  (case pre-kind
                        [(0) (and same-cc (fx<=? distance2 400))]  ;; p6l4 axon span at l4
                        [(1) (fx<=? distance2 9)]                  ;; p23 axon span at l4
                        [(2) (and same-cc (fx<=? distance2 500))]  ;; ss4l4 axons intra-cc only
                        [(3) (and same-cc (fx<=? distance2 100))]  ;; ss4l23 axon span at l4
                        [(4) (fx<=? distance2 9)]) ]               ;; p4 axon span
                    [(p4a)  (case pre-kind
                        [(0) (fx<=? distance2 9)]                  ;; p6l4 axon span at l2
                        [(1) (fx<=? distance2 500)]                ;; p23 axon span at l2
                        [(2) (and #;same-cc (fx<=? distance2 36))]   ;; ss4l4 axon span at l2
                        [(3) (and same-cc (fx<=? distance2 64))]   ;; ss4l23 axon span at l2
                        [(4) (fx<=? distance2 500)])]              ;; p4 axon span at l2
                    )])
              (define (connect? ss4l4 ss4l23 p4 apical)
                ;; randomise connect? by proportions ss4l4:ss4l23:p4
                (let* ( 
                    [split1 (fx* ss4l4 ss4l4-c/mc)]
                    [split2 (fx+ split1 (fx* ss4l23 ss4l23-c/mc))]
                    [target (random (fx+ split2 (fx* p4 p4-c/mc)))])
                  (case post-kind
                    [(ss4l4)   (fx<? target split1)]
                    [(ss4l23)  (fx<? (fx1- split1) target split2)]
                    [(p4b)     (fx<=? split2 target)]
                    [(p4a)     apical] )))
              ;; Modified Peters' rule: p6l4 axons prefer ss4l23+p4, ss4l4 axons prefer ss4l4
              (case pre-kind
                [(0)    (and reach? (connect? 25 75 75 #f)) ]         ;; 25 75 75
                [(1)    reach? ]
                [(2)    (and reach? (connect? 75 25 25 #f)) ]         ;; 75 25 25
                [(3 4)  (and reach? (connect? 25 75 75 #t)) ] ))))))  ;; 25 75 75
                                                                                            ;
    (assert (fx<=? ncc (expt 2 ccx-bits)))
    (assert (fx<=? (fxmax
        (fx* nmc (fxmax ss4l4-c/mc ss4l23-c/mc p4-c/mc l6-c/mc l2-c/mc)))
        (expt 2 cellx-bits)))
    (new
      (build-vector ncc (lambda (ccx)
          (l2:make-l2 ncc ccx nmc l2-c/mc cp-overrides l2-connect?)))
      (build-vector ncc (lambda (ccx)
          (l4:make-l4 ccx nmc ss4l4-c/mc ss4l23-c/mc p4-c/mc attm-overrides l4-connect?)))
      #;(build-vector ncc (lambda (ccx)
          (if (fx=? nmc minicolumns/macrocolumn)             ;; descending sort so that sources-for-ccs sorts up
            (sort! fx>? (filter (lambda (ccy)
                (fx<? (cc-distance2 ccx ccy) 800))
              (iota ncc)))
            (list ccx))))                ;; just the cc if not using hexagonal topology
      (build-vector ncc (lambda (ccx)
          (if (fx=? nmc minicolumns/macrocolumn)             ;; descending sort so that sources-for-ccs sorts up
            (sort! fx>?
              (if (fx>=? ncc 7)          ;; cc and 3 others for patchy connectivity
                (cons ccx (u32-sample (remv ccx (iota ncc)) 3))
                (iota ncc)))
            (list ccx))))                ;; just the cc if not using hexagonal topology
      )))))
                                                                                            ;
(define (compute p features locations    ;; Patch (CCVecOf SDR) (CCVecOf SDR) Boolean ->
          learn)
  ;; run one timestep of patch
  (define (sources ccx layer cellxs)     ;; Nat Nat {CellX} -> {Source}
    ;; produce source list for one cc/layer; creates new pairs so not changed by compute
    (let ([stride (make-source ccx layer 0)])
      (map (lambda (cellx) (fx+ cellx stride))
        cellxs)))
  (define (all-sources proc)             ;; (CCX L4 -> { [Nat . {CellX}] }) -> {Source}
    ;; produce sorted list of sources for all ccs and specified layers/cells
    (do ( [ccx (fx1- (vector-length (patch-l4s p))) (fx1- ccx)]
          [acc (list) (append!
                        (fold-right (lambda (s acc)
                            (append! (sources ccx (car s) (cdr s)) acc))
                          (list)
                          (proc ccx (vector-ref (patch-l4s p) ccx)))
                        acc)] )
        ((fxnegative? ccx) acc) ) )
  (define (sources-for-ccs proc)         ;; (CCX L4 -> { [Nat . {CellX}] }) -> (CCVecOf {Source})
    ;; like all-sources, but produce cc vector for nearby ccs and specified layers/cells
    (vector-map (lambda (ccx)
        (do ( [ccxx (vector-ref (patch-nearby-ccs p) ccx) (cdr ccxx)]
              [acc (list) (append!
                            (fold-right (lambda (s acc)
                                (append! (sources (car ccxx) (car s) (cdr s)) acc))
                              (list)
                              (proc (car ccxx) (vector-ref (patch-l4s p) (car ccxx))))
                            acc)] )
            ((null? ccxx) acc) ) )
      (indexes (patch-nearby-ccs p))))
  ;; (compute)
  ;; assemble all L4 inputs (from nearby ccs for each cc) before cc parallel depolarize-cells
  ;; inputs to depolarize are active cells from previous timestep
  (let (
      [locations
        (sources-for-ccs (lambda (ccx l4)
            `([0 . ,(vector-ref locations ccx)]))) ]
      [l4-activitys               ;; (CCVecOf {Source})
        (sources-for-ccs (lambda (ccx l4)
            `([2 . ,(l4:get-active-cells l4 'ss4l4)]
              [3 . ,(l4:get-active-cells l4 'ss4l23)]
              [4 . ,(l4:get-active-cells l4 'p4)]))) ]
      [l2-activitys
        (sources-for-ccs (lambda (ccx l4)
            `([1 . ,(l2:get-active-cells (vector-ref (patch-l2s p) ccx))]
              [3 . ,(l4:get-active-cells l4 'ss4l23)]
              [4 . ,(l4:get-active-cells l4 'p4)]))) ] )
    (threaded-vector-for-each            ;; thread per cortical column
      (lambda (l4 location l4-activity l2-activity)
        (l4:depolarize-cells l4 
          location
          l4-activity
          l2-activity
          learn))
      (patch-l4s p) locations l4-activitys l2-activitys)
    ;; depolarized L4 cells with tc input fire first: projection to other cols could inhibit?
    (vector-for-each (lambda (l4 feature location)
        (l4:inhibit-columns l4 feature location))
      (patch-l4s p) features locations)
    (let (
        [l4-learnings
          (sources-for-ccs (lambda (ccx l4)
              `([2 . ,(l4:get-learning-cells l4 'ss4l4)]
                [3 . ,(l4:get-learning-cells l4 'ss4l23)]
                [4 . ,(l4:get-learning-cells l4 'p4)]))) ]
        [l2-learnings
          (sources-for-ccs (lambda (ccx l4)
              `([1 . ,(l2:get-active-cells (vector-ref (patch-l2s p) ccx))]
                [3 . ,(l4:get-learning-cells l4 'ss4l23)]
                [4 . ,(l4:get-learning-cells l4 'p4)]))) ] )
      ;; with inhibited cols/learning cells for all ccs known, activate-cells in parallel    
      (threaded-vector-for-each            ;; thread per cortical column
        (lambda (l4 feature location l4-activity l4-learning l2-activity l2-learning)
          (l4:activate-cells l4 
            feature                        ;; proximal input
            location
            l4-activity
            l4-learning
            l2-activity
            l2-learning
            learn))
        (patch-l4s p) features locations  l4-activitys l4-learnings l2-activitys l2-learnings)))
  ;; L2 uses inhibited minicols from L4; with hexagonal topology distal input includes laterals
  (let (
      [l4-activitys
        (sources-for-ccs (lambda (ccx l4)
            `([3 . ,(l4:get-active-cells l4 'ss4l23)]
              [4 . ,(l4:get-active-cells l4 'p4)]))) ]
      [l2-activitys
        (sources-for-ccs (lambda (ccx l4)
            `([1 . ,(l2:get-active-cells (vector-ref (patch-l2s p) ccx))])))]
      [l4-predicteds
        (sources-for-ccs (lambda (ccx l4)   
            `([3 . ,(l4:get-predicted-cells l4 'ss4l23)]
              [4 . ,(l4:get-predicted-cells l4 'p4)]))) ]
      [l2-laterals                     ;; lateral inputs not used if nearby ccs 
        (if (fx>? (length (vector-ref (patch-nearby-ccs p) 0)) 1)
          (make-vector (vector-length (patch-nearby-ccs p)) '())
          (let ([l2s  (patch-l2s p)])
            (vector-map (lambda (ccx)
                (vector-fold-left (lambda (acc ccy)
                    (if (fx=? ccy ccx)  acc
                      (cons (sources ccy 1 (l2:get-active-cells (vector-ref l2s ccy))) acc)))
                (list)
                (indexes l2s)))
              (indexes l2s)))) ])
    (threaded-vector-for-each            ;; thread per cortical column
      (lambda (l2 l4 l4-activity l2-activity l2-lateral l4-predicted)
        (l2:compute l2
          l4-activity                    ;; feedforward-input
          l2-activity                    ;; distal-input
          l2-lateral                     ;; lateral-inputs
          l4-activity                    ;; feedforward-growth-candidates
          learn
          l4-predicted                   ;; predicted-input
          (l4:l4-bursting-cols l4)))     ;; inhibited-cols
      (patch-l2s p) (patch-l4s p) l4-activitys l2-activitys l2-laterals l4-predicteds)))
                                                                                            ;
(define (reset p)
  (vector-for-each  l2:reset  (patch-l2s p))
  (vector-for-each  l4:reset  (patch-l4s p)))
                                                                                            ;
(define (reset-l2 p)
  (vector-for-each  l2:reset  (patch-l2s p)))
                                                                                            ;
(define (reset-l4 p)
  (vector-for-each  l4:reset  (patch-l4s p)))
                                                                                            ;
(define (reset-seq p)
  (vector-for-each  l4:reset-seq  (patch-l4s p)))
                                                                                            ;
(define (print-statistics p)
  (l2:print-statistics
    (lambda (f l)
      (vector-fold-left (lambda (sum l2)
          (+ sum (f (l l2))))
        0
        (patch-l2s p))))
  (l4:print-statistics (patch-l4s p)))

)