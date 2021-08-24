#lang racket

(require plot)

(define (with-x-coords vars)            ;; (listof (listof Number)) -> (listof (listof (vector x y))
  ;; each element of vars is a list of y coordinates for a plot line
  (map (lambda (y-coords)
         (if y-coords
             (map vector (range (length y-coords)) y-coords)
             '()))
       vars))

(define (evens vs)         ;; (listof (vector x y))
  (filter (lambda (v) (even? (vector-ref v 0)))
          vs))
   
(define (odds vs)         ;; (listof (vector x y))
  (filter (lambda (v) (odd? (vector-ref v 0)))
          vs))

(define (up vs)
  (map (lambda (v)
         (vector (vector-ref v 0)
                 (if (odd? (vector-ref v 0))
                     (+ (vector-ref v 1) 0.25)
                     (- (vector-ref v 1) 0.25))))
       vs))
  
(define (down vs)
  (map (lambda (v)
         (vector (vector-ref v 0)
                 (if (even? (vector-ref v 0))
                     (+ (vector-ref v 1) 0.5)
                     (- (vector-ref v 1) 0.5))))
       vs))
  
(define (plus-lines vs #:color c #:alpha [a 1.0] . params)
  (list
   (keyword-apply
    lines (car params) (cadr params)
    vs '() #:color c #:alpha a)
   (keyword-apply
    points '(#:line-width #:size #:sym) '(2.5 8 plus)
    vs '() #:color c #:alpha a)))

(define x-max (make-parameter 0))
(define y-max (make-parameter 0))

(define (offset vs ws by)    ;; (listof (vector x y))
  (if (= (vector-ref (car vs) 1) (vector-ref (car ws) 1))
      (map (lambda (v)
             (vector (+ (* .0025 by (x-max)) (vector-ref v 0)) (+ (* .005 by (y-max)) (vector-ref v 1))))
           vs)
      vs))

(define render
  (case-lambda
    [(dp4 pp4 dss4l4 pss4l4) (render dp4 pp4 dss4l4 pss4l4 '() '()) ]
    [(dp4 pp4 dss4l4 pss4l4 dss4l23 pss4l23)
     (append
      (list
       (plus-lines (offset dp4 pp4 1) #:color 'RoyalBlue #:alpha 0.5
                   '(#:label #:style) '("Depolarized sensorimotor (p4)" short-dash))
       (plus-lines (offset pp4 pss4l23 1) #:color 'RoyalBlue
                   '(#:label) '("Correctly predicted sensorimotor (p4)"))
       (if (pair? dss4l23)
           (list
            (plus-lines (offset dss4l23 pss4l23 1) #:color 'Red #:alpha 0.5
                        '(#:label #:style) '("Depolarized sensorimotor ss4(L2/3)" short-dash))
            (plus-lines (offset pss4l23 pp4 (- 1)) #:color 'Red
                        '(#:label) '("Correctly predicted sensorimotor ss4(L2/3)")))
           '())
       (plus-lines (offset dss4l4 pss4l4 1) #:color 'MediumSeaGreen #:alpha 0.5
                   '(#:label #:style) '("Depolarized sequence ss4(L4)" short-dash))
       (plus-lines (offset pss4l4 dss4l4 (- 1)) #:color 'MediumSeaGreen
                   '(#:label) '("Correctly predicted sequence ss4(L4)")))) ]))

(define render6
  (case-lambda
    [(psm pts) (render6 psm pts '()) ] 
    [(psm pts ptx)
     (append
      (if (pair? ptx)
          (list (lines (offset ptx psm (- 0.5)) #:color 'DarkOrange
                       #:label "Predicted active sensorimotor ss4(L2/3) cells"))
          '())
      (list (lines (offset psm ptx 0.5) #:color 'Red
                   #:label "Predicted active sensorimotor p4 cells"))
      (list (lines pts #:color 'RoyalBlue
                   #:label "Predicted active temporal seq ss4(L4) cells"))) ]))

(define (renderH3 rs ys xs)
  (let* ((final   (for/list ([y ys] [x xs] #:when   (>= x 8)) y))
         (correct (for/list ([f final]     #:when   (member f rs)) f))
         (missed  (append
                   (for/list ([f final]    #:unless (member f correct)) f)
                   (for/list ([r rs]       #:unless (member r correct)) r))))
    (list
     (map (lambda (y x)
            (hrule y 0 x #:width 1 #:style 'solid))
          ys xs)
     (map (lambda (x)
            (vrule x 30 4090 #:width 1 #:color 'white))
          (range 1 15))
     (points (map vector (make-list (length correct) 15)
                  (map (lambda (c) (+ c 4)) correct))) ;; centre circle on line
     (points (map vector (make-list (length missed) 15) missed) #:sym 'times))))
         
(define (clear-legend ys)
  (+ (* 5 (length ys))
     (apply max
            (apply append
                   (map (lambda (y) (if y (cddr y) '(0))) ys)))))

(define plots '())

(define (plot-experiment data)
  (define (value-for key)
    (let ((entry (assoc key data)))
      (if entry (if (pair? (cdr entry)) (cadr entry) (cdr entry)) #f)))
  (let ( (figure (value-for "figure"))
         (using  (value-for "using"))
         (l2r    (value-for "l2r"))
         (l2r1   (value-for "l2r1"))
         (l2r2   (value-for "l2r2"))
         (l2a    (value-for "l2a"))
         (l2ac   (value-for "l2ac"))
         (l2a1   (value-for "l2a1"))
         (l2a1c  (value-for "l2a1c"))
         (l2a2   (value-for "l2a2"))
         (l2a2c  (value-for "l2a2c"))
         (l4lnp  (value-for "l4lnp"))
         (l3pa   (value-for "l3pa"))
         (l3lpa  (value-for "l3lpa"))
         (l4pa   (value-for "l4pa"))
         (l4lpa  (value-for "l4lpa"))
         (TMlnp  (value-for "TMlnp"))
         (TMpa   (value-for "TMpa"))
         (TMlpa  (value-for "TMlpa"))
         (TXlnp  (value-for "TXlnp"))
         (TXpa   (value-for "TXpa"))
         (TXlpa  (value-for "TXlpa")))
    (when (and using (null? plots))
      (for-each
       (lambda (u)
         (display " '")
         (display u))
       using))
    (parameterize
        ([plot-width  450]
         [plot-height 350]
         [plot-font-size 12]
         [plot-font-face "Helvetica"]
         [plot-legend-anchor 'top-right]
         [plot-x-ticks (linear-ticks #:number 10 #:divisors '(1))]
         [plot-x-far-ticks no-ticks]
         [plot-y-far-ticks no-ticks]
         [plot-x-label #f] [plot-y-label #f]
         [line-width 3])
      (set!
       plots
       (append
        plots
        (list
         (case figure
           [("kt1")
            (let ((course (value-for "course")))
              (plot (lines #:width 1 #:x-min 0 #:x-max 256 #:y-min 0 #:y-max 256
                           (map vector
                                (map real-part course)
                                (map imag-part course)))
                    #:width 512 #:height 512))]
           [("kt2" "kt3" "kt4")
            (let ((enc-rf (value-for "enc-rf")))
              (parameterize ([plot-x-tick-labels? #f]
                             [plot-y-tick-labels? #f]
                             [plot-tick-size 0])
                (plot (points enc-rf #:sym 'dot
                              #:x-min 0 #:x-max 255 #:y-min 0 #:y-max 255)
                      #:width 130 #:height 130)))]
           [("f4a")
            (let ((ys (if TXlnp
                          (list l4lnp l4pa TMlnp TMpa TXlnp TXpa)
                          (list l4lnp l4pa TMlnp TMpa))))
              (parameterize ([plot-y-ticks (linear-ticks #:number 6 #:divisors '(5))]
                             [x-max 9.5] [y-max (clear-legend ys)])
                (plot
                 #:title "Figure 4A' Average of predictions inferring 50 sequences"
                 #:x-min -0.5 #:x-max (x-max) #:x-label "Input number"
                 #:y-min -5   #:y-max (y-max) #:y-label "Number of cells"
                 (apply render (with-x-coords ys)))))]
           [("f5a")
            (let ((ys (if TXlnp
                          (list l4lnp l4pa TMlnp TMpa TXlnp TXpa)
                          (list l4lnp l4pa TMlnp TMpa))))
              (parameterize ([plot-y-ticks (linear-ticks #:number 6 #:divisors '(5))]
                             [x-max 9.5] [y-max (clear-legend ys)])
                (plot 
                 #:title "Figure 5A' Average of predictions inferring 50 objects"
                 #:x-min -0.5 #:x-max (x-max) #:x-label "Input number"
                 #:y-min -5   #:y-max (y-max) #:y-label "Number of cells"
                 (apply render (with-x-coords ys)))))]

           [("f6")
            (let ([so-y (+ 1 (* 1.1 (apply max (append
                                                (apply append
                                                       (if l4pa (map vector->list l4pa) '()))
                                                (apply append
                                                       (if TMpa (map vector->list TMpa) '()))
                                                (apply append
                                                       (if TXpa (map vector->list TXpa) '()))))))])
              (define (transpose vs)   ;; (listof (vectorof Number)) -> (listof (listof Number))
                ;; all first elements of vectors as list, then all second, etc
                (build-list
                 (vector-length (car vs))
                 (lambda (i)
                   (map (lambda (v)
                          (if (negative? (vector-ref v i))
                              (- so-y 1)
                              (vector-ref v i))
                          #;(vector-ref v i))
                        vs))))
              (define (with-x-coords y-coords)            ;; (listof Number) -> (listof (vector x y))
                (map vector (range (length y-coords)) y-coords))
              (define (f6lines lpa1 lpa2 os col lab)
                (if (vector? (car lpa1))
                    (list
                     (lines '(#(0 0) #(0 0)) #:color col #:alpha 0.8 #:label lab)
                     (map (lambda (l1 l2)
                            (lines (offset (with-x-coords l1) (with-x-coords l2) os)
                                   #:color col #:alpha (max 0.25 (/ 0.8 (vector-length (car lpa1))))))
                          (transpose lpa1)(transpose lpa2)))
                    (lines (with-x-coords lpa1) #:color col  #:label lab)))
            
              (parameterize
                  ([plot-width  700]
                   [plot-height 450]
                   [plot-font-size 14]
                   [plot-x-ticks (linear-ticks #:number 21 #:divisors '(2))]
                   [plot-y-ticks (linear-ticks #:number 6 #:divisors '(1))]
                   [x-max 101]
                   [y-max (+ so-y (/ so-y 3))])
                (plot 
                 #:title "Figure 6a Combined sensorimotor/sequence streams"
                 #:x-min -2 #:x-max (x-max) #:x-label "Input number"
                 #:y-min -1 #:y-max (y-max) #:y-label "Number of cells"
                 #;(list  ;; plot number predicted
                    (if (pair? TXlpa)
                        (f6lines TXlpa l4lpa -0.5 'DarkOrange "Predicted active sensorimotor ss4(L2/3) cells")
                        '())
                    (if (pair? l3lpa)
                        (f6lines l3lpa l4pa -0.5 'Green "Predicted active temporal seq p3 cells")
                        '())
                    (f6lines l4lpa TXlpa 0.5 'Red "Predicted active sensorimotor p4 cells")
                    (f6lines TMlpa TMlpa 0.0 'RoyalBlue "Predicted active temporal seq ss4(L4) cells")
                    (map (lambda (x)
                           (vrule x 0 so-y #:width 1 #:style 'long-dash))
                         (range -0.5 90 10))
                    (map (lambda (x l)
                           (point-label (vector x so-y) l #:size 12 #:point-size 0))
                         (range 0 100 10)
                         (let ((s "Sequence") (o "   Object"))
                           (list s o s o s s o s o o))))
                 (list  ;; plot correctly predicted
                  (if (pair? TXpa)
                      (f6lines TXpa (if (pair? l4pa) l4pa TXpa)
                               -0.5 'DarkOrange "Predicted active sensorimotor ss4(L2/3) cells")
                      '())
                  (if (pair? l3pa)
                      (f6lines l3pa TXpa 0.5 'Green "Predicted active temporal seq p3 cells")
                      '())
                  (if (pair? l4pa)
                      (f6lines l4pa (if (pair? TXpa) TXpa l4pa)
                               0.5 'Red "Predicted active sensorimotor p4 cells")
                      '())
                  (if (pair? TMpa)
                      (f6lines TMpa TMpa 0.0 'RoyalBlue "Predicted active temporal seq ss4(L4) cells")
                      '())
                  (map (lambda (x)
                         (vrule x 0 so-y #:width 1 #:style 'long-dash))
                       (range -0.5 90 10))
                  (map (lambda (x l)
                         (point-label (vector x so-y) l #:size 12 #:point-size 0))
                       (range 0 100 10)
                       (let ((s "Sequence") (o "   Object"))
                         (list s o s o s s o s o o))))
                 )))]

           [("H3b" "H3c" "H4b")
            (parameterize
                ([plot-width  200]
                 [plot-height 500]
                 [plot-x-ticks (linear-ticks #:number 3 #:divisors '(2))]
                 [plot-y-ticks (linear-ticks #:number 9 #:divisors '(2))]
                 [plot-tick-size 0]
                 [plot-x-far-axis? #f]
                 [plot-y-far-axis? #f])
              (let ((one-col (string=? figure "H3b")))
                (display (plot
                          #:title "   Column 1"
                          #:x-min 0 #:x-max 15.3 #:x-label (if one-col "Number of sensations" "")
                          #:y-min -20 #:y-max 4116 #:y-label "Neuron #"
                          (renderH3 l2r l2a l2ac)))
                (unless one-col
                  (parameterize ([plot-y-axis? #f] [plot-width 150])
                    (display (plot
                              #:title "   Column 2"
                              #:x-min 0 #:x-max 15.3 #:x-label "Number of sensations"
                              #:y-min -20 #:y-max 4116 #:y-label #f
                              (renderH3 l2r1 l2a1 l2a1c)))
                    (display (plot
                              #:title "   Column 3"
                              #:x-min 0 #:x-max 15.3 #:x-label ""
                              #:y-min -20 #:y-max 4116 #:y-label #f
                              (renderH3 l2r2 l2a2 l2a2c)))))
                (newline)
                (for-each
                 (lambda (ac r)
                   (display (length r)) (display ": ")
                   (display
                    (map
                     (lambda (n)
                       (length
                        (filter (lambda (x) (> x n)) ac)))
                     (range 20)))
                   (newline))
                 (if one-col (list l2ac) (list l2ac l2a1c l2a2c))
                 (if one-col (list l2r) (list l2r l2r1 l2r2)))
                " "))]
           ))))
      ;(for-each display `("\n" ,figure ": "))
      )))

(define (newest-data)
  (argmax file-or-directory-modify-seconds
          (for/fold
           ([data-files (list)])
           ([d (directory-list)])
            (if (directory-exists? d)
                (let ((f (build-path d "experiment.data")))
                  (if (file-exists? f)
                      (cons f data-files)
                      data-files))
                data-files))))

(let loop ((modified 0))
  (let* ((f (newest-data))
         (t (file-or-directory-modify-seconds f)))
    (when (> t modified)
      (with-input-from-file f
        (lambda ()
          (let loop ([data (read)])
            (file-position (current-input-port) 0)
            (when (eof-object? data)
              (sleep 1)
              (loop (read)))
            (newline)
            (display f)
            (set! plots '())
            (let plot ([data data])
              (plot-experiment (car data))
              (unless (null? (cdr data)) (plot (cdr data))))
            (newline)
            (if #f #;(> (length plots) 5) (begin
                                            (display (car plots))
                                            (newline)
                                            (for-each display (cdr plots)))
                (for-each display plots))
            (newline)
            (sleep 1)
            (sync (filesystem-change-evt f))
            (loop (read))))))
    (loop t)))
