#lang racket

(require plot)

(define (with-x-coords vars)            ;; (listof (listof Number))
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

(define (render dsm psm dts pts)
  (list
   (plus-lines (offset dsm psm 1) #:color 'Red #:alpha 0.5
               '(#:label #:style) '("Depolarized sensorimotor cells" short-dash))
   (plus-lines (offset psm dsm (- 1)) #:color 'RoyalBlue
               '(#:label) '("Predicted active sensorimotor cells"))
   (plus-lines (offset dts pts 1) #:color 'DarkOrange #:alpha 0.5
               '(#:label #:style) '("Depolarized temporal sequence cells" short-dash))
   (plus-lines (offset pts dts (- 1)) #:color 'MediumSeaGreen
               '(#:label) '("Predicted active temporal sequence cells"))))

(define (render6 psm pts)
  (list
   (lines psm #:color 'Red
          #:label "Predicted active cells in sensorimotor layer")
   (lines pts #:color 'RoyalBlue
          #:label "Predicted active cells in temporal sequence layer")))

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
  (+ (* 10 (length ys))
     (apply max
            (apply append
                   (map cddr ys)))))
  
(define (plot-experiment data)
  (define (value-for key)
    (let ((entry (assoc key data)))
      (if entry  (cadr entry)  #f)))
  (let ( (figure (value-for "figure"))
         (using  (value-for "using"))
         (L2r    (value-for "L2r"))
         (L2r1   (value-for "L2r1"))
         (L2r2   (value-for "L2r2"))
         (L2a    (value-for "L2a"))
         (L2ac   (value-for "L2ac"))
         (L2a1   (value-for "L2a1"))
         (L2a1c  (value-for "L2a1c"))
         (L2a2   (value-for "L2a2"))
         (L2a2c  (value-for "L2a2c"))
         (L4lnp  (value-for "L4lnp"))
         (L4pa   (value-for "L4pa"))
         (L4lpa  (value-for "L4lpa"))
         (TMlnp  (value-for "TMlnp"))
         (TMpa   (value-for "TMpa"))
         (TMlpa  (value-for "TMlpa")))
    (parameterize
        ([plot-width  450]
         [plot-height 350]
         [plot-font-size 12]
         [plot-font-face "Helvetica"]
         [plot-legend-anchor 'top-right]
         [plot-x-ticks (linear-ticks #:number 10 #:divisors '(1))]
         [plot-x-far-ticks no-ticks]
         [plot-y-far-ticks no-ticks]
         [line-width 3])
      (display
       (case figure
         [("f4a")
          (let ((ys (list L4lnp L4lpa TMlnp TMlpa)))
            (parameterize ([plot-y-ticks (linear-ticks #:number 6 #:divisors '(5))]
                           [x-max 9.5] [y-max (clear-legend ys)])
              (plot
               #:title "Figure 4A' Average of predictions inferring 50 sequences"
               #:x-min -0.5 #:x-max (x-max) #:x-label "Input number"
               #:y-min -5   #:y-max (y-max) #:y-label "Number of cells"
               (apply render (with-x-coords ys)))))]
         [("f5a")
          (let ((ys (list L4lnp L4lpa TMlnp TMlpa)))
            (parameterize ([plot-y-ticks (linear-ticks #:number 6 #:divisors '(5))]
                           [x-max 9.5] [y-max (clear-legend ys)])
              (plot 
               #:title "Figure 5A' Average of predictions inferring 50 objects"
               #:x-min -0.5 #:x-max (x-max) #:x-label "Input number"
               #:y-min -5   #:y-max (y-max) #:y-label "Number of cells"
               (apply render (with-x-coords ys)))))]
         [("f6")
          (parameterize
              ([plot-width  600]
               [plot-height 450]
               [plot-font-size 14]
               [plot-x-ticks (linear-ticks #:number 17 #:divisors '(2))]
               [plot-y-ticks (linear-ticks #:number 6 #:divisors '(1))])
            (let* ((ys (list L4lpa TMlpa))
                   (so-y (apply max (apply append ys)))
                   (so-y (+ so-y (/ so-y 20))))
              
              (plot 
               #:title #;"Figure 6' Inferring combined sensorimotor and temporal sequence stream"
               "Figure 6' Combined sensorimotor/sequence streams, 7 cortical columns of 150 minicolumns"
               #:x-min -2 #:x-max 81 #:x-label "Input number"
               #:y-min -1 #:y-max (+ so-y (/ so-y 4)) #:y-label "Number of cells"
               (list
                (apply render6 (with-x-coords ys))
                (map (lambda (x)
                       (vrule x 0 (+ 0 so-y) #:width 1 #:style 'long-dash))
                     (range -0.5 70 10))
                (map (lambda (x l)
                       (point-label (vector x so-y) l #:size 12 #:point-size 0))
                     (range 0 80 10)
                     (let ((s "Sequence") (o "   Object"))
                       (list s o s o s s o s)))))))]
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
                        (renderH3 L2r L2a L2ac)))
              (unless one-col
                (parameterize ([plot-y-axis? #f] [plot-width 150])
                  (display (plot
                            #:title "   Column 2"
                            #:x-min 0 #:x-max 15.3 #:x-label "Number of sensations"
                            #:y-min -20 #:y-max 4116 #:y-label #f
                            (renderH3 L2r1 L2a1 L2a1c)))
                  (display (plot
                            #:title "   Column 3"
                            #:x-min 0 #:x-max 15.3 #:x-label ""
                            #:y-min -20 #:y-max 4116 #:y-label #f
                            (renderH3 L2r2 L2a2 L2a2c)))))
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
               (if one-col (list L2ac) (list L2ac L2a1c L2a2c))
               (if one-col (list L2r) (list L2r L2r1 L2r2)))
              " "))]
         ))
      (display "\n(experiment-")
      (display figure)
      (for-each
       (lambda (u)
         (display " '")
         (display u))
       using)
      (display ")\n")
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
          (display f) (newline)
          (plot-experiment (read)))))
    (sleep 2)
    (loop t)))







#;(with-input-from-file f
    (lambda ()
      (let loop ((data (read)))
        (file-position (current-input-port) 0)
        (if (eof-object? data) (exit)
            (begin
              (plot-cs data)
              (let wait ()
                (sleep 0.1)
                (sync (filesystem-change-evt f))
                (let ((newdata (read)))
                  (file-position (current-input-port) 0)
                  (if (equal? newdata data)
                      (wait)
                      (loop newdata)))))))))
