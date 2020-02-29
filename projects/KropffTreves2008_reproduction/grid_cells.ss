#!chezscheme

(library (HTM-scheme projects KropffTreves2008_reproduction grid_cells)

(export
  make-gc
  reset
  compute)
  
(import
  (except (chezscheme) add1 make-list reset)
  (except (HTM-scheme HTM-scheme algorithms htm_prelude) random)
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
  
(define (compute gc input-vector learn)  ;; GC InputVec Boolean -> (listof ColumnX)
  ;; produce active columns from input; optionally update sp if learning
  (sp:sp-iteration-num-set! gc (fx1+ (sp:sp-iteration-num gc)))
  (when learn (sp:sp-iteration-learn-num-set! gc (fx1+ (sp:sp-iteration-learn-num gc))))
  (let* (
      [overlaps         (sp:calculate-overlap gc input-vector)]  ;; (ColVecOf Overlap)
      [sum-synapses     (sp:sp-connected-counts gc)]             ;; (ColVecOf Nat) set by adapt-synapses
      [n-overlaps       (vector-map (lambda (ov ss)
                            (if (fxzero? ss)  0.0 
                                (fl/ ov (fixnum->flonum ss))))
                          overlaps sum-synapses)]
      [f-overlaps       (_fatigue gc n-overlaps)]
      [boosted-overlaps (if learn
                            (vector-map fl* f-overlaps (sp:sp-boost-factors gc))
                            n-overlaps)]
      [active-columns   (sp:inhibit-columns gc boosted-overlaps)])

    (when learn
      (adapt-synapses          gc input-vector active-columns)
      (sp:update-duty-cycles   gc overlaps active-columns)
      (sp:bump-up-weak-columns gc)
      (sp:update-boost-factors gc)
      (when (fxzero? (fxmod (sp:sp-iteration-num gc) (sp:sp-update-period gc)))
          (sp:sp-inhibition-radius-set! gc (sp:update-inhibition-radius gc))
          (sp:update-min-duty-cycles gc)))
    active-columns))
    
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

(define (adapt-synapses                  ;; SP InputVec (listof ColumnX) ->
          gc input-vector active-columns)
  ;; update permanences in segments of active columns (+ if synapse's input on, - if not)
  (let ((syn-perm-trim-threshold (sp:sp-syn-perm-trim-threshold gc))
        (syn-perm-active-inc     (sp:sp-syn-perm-active-inc gc))
        (syn-perm-inactive-dec   (sp:sp-syn-perm-inactive-dec gc))
        (syn-perm-connected      (sp:sp-syn-perm-connected gc)))
    (sp:for-each-segment gc (lambda (segment cx)
        (let ((prev-active (memv cx (gc-prev-active-columns gc))))
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
                          (synapses-set! segment i
                            (sp:increase-perm synapse syn-perm-trim-threshold)))]
                      [(fx=? perm max-perm)
                        (unless inp-bit
                          (synapses-set! segment i
                            (sp:decrease-perm synapse syn-perm-inactive-dec syn-perm-trim-threshold)))]
                      #;[else (synapses-set! segment i
                              (if inp-bit
                                (sp:increase-perm synapse syn-perm-active-inc)
                                (sp:decrease-perm synapse syn-perm-inactive-dec syn-perm-trim-threshold)))])))))
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
          ))
      active-columns))
  (gc-prev-active-columns-set! gc active-columns)
  (gc-prev-input-vector-set!   gc input-vector))
                                                                              ;
)
