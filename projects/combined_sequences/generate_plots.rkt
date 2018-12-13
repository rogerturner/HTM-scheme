#lang racket

(require plot)

(define (plus-lines vs #:color c #:label l)
  (list
   (lines vs #:color c #:label l)
   (keyword-apply
    points '(#:line-width #:size #:sym) '(2.5 8 plus)
    vs '() #:color c)))

(define (render l4pa tmnp tmpa)
  (list
   (plus-lines l4pa #:color 'RoyalBlue
               #:label "Predicted active cells in sensorimotor layer")
   (plus-lines tmnp #:color 'DarkOrange
               #:label "Predicted cells in temporal sequence layer")
   (plus-lines tmpa #:color 'MediumSeaGreen
               #:label "Predicted active cells in temporal sequence layer")))

(define (render6 l4pa tmpa)
  (list
   (lines l4pa #:color 'Red
          #:label "Predicted active cells in sensorimotor layer")
   (lines tmpa #:color 'RoyalBlue
          #:label "Predicted active cells in temporal sequence layer")))

(define (renderH3 rs ys xs)
  (list
   (points (map vector (make-list (length rs) 15)
                (map (lambda (r) (+ r 4)) rs))) ;; tweak to match line centre
   (map (lambda (y x)
          (hrule y 0 x #:width 2 #:style 'solid))
        ys xs)
   (map (lambda (x)
          (vrule x 30 4090 #:width 1 #:color 'white))
        (range 1 15))))
         
(define (with-x-coords n vars)            ;; Nat (listof (listof Number))
  ;; each element of vars is a list of y coordinates for a plot line
  (map (lambda (y-coords)
         (if y-coords
             (map vector (range n) y-coords)
             '()))
       vars))

(define (clear-legend ys)
  (+ (* 10 (length ys))
     (apply max
            (apply append
                   (map cddr ys)))))
  
(define (plot-cs data)
  (define (value-for key)
    (let ((entry (assoc key data)))
      (if entry  (cadr entry)  #f)))
  (let ( (figure (value-for "figure"))
         (L2r    (value-for "L2r"))
         (L2r1   (value-for "L2r1"))
         (L2r2   (value-for "L2r2"))
         (L2a    (value-for "L2a"))
         (L2ac   (value-for "L2ac"))
         (L2a1   (value-for "L2a1"))
         (L2a1c  (value-for "L2a1c"))
         (L2a2   (value-for "L2a2"))
         (L2a2c  (value-for "L2a2c"))
         (L4lpa  (value-for "L4lpa"))
         (TMlnp  (value-for "TMlnp"))
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
      (case figure
        [("A4a")
         (parameterize ([plot-y-ticks (linear-ticks #:number 6 #:divisors '(5))])
           (let ((ys (list L4lpa TMlnp TMlpa)))
             (plot
              #:title "Figure 4A' Average of predictions inferring 50 sequences"
              #:x-min -0.5 #:x-max 9.5 #:x-label "Input number"
              #:y-min -5   #:y-max (clear-legend ys) #:y-label "Number of cells"
              (apply render (with-x-coords 10 ys)))))]
        [("A5a")
         (parameterize ([plot-y-ticks (linear-ticks #:number 6 #:divisors '(5))])
           (let ((ys (list L4lpa TMlnp TMlpa)))
             (plot 
              #:title "Figure 5A' Average of predictions inferring 50 objects"
              #:x-min -0.5 #:x-max 9.5 #:x-label "Input number"
              #:y-min -5   #:y-max (clear-legend ys) #:y-label "Number of cells"
              (apply render (with-x-coords 10 ys)))))]
        [("A6")
         (parameterize
             ([plot-width  700]
              [plot-height 500]
              [plot-font-size 16]
              [plot-x-ticks (linear-ticks #:number 17 #:divisors '(2))]
              [plot-y-ticks (linear-ticks #:number 6 #:divisors '(1))])
           (let ((ys (list L4lpa TMlpa)))
             (plot 
              #:title "Figure 6' Inferring combined sensorimotor and temporal sequence stream"
              #:x-min -2 #:x-max 81 #:x-label "Input number"
              #:y-min -1 #:y-max (clear-legend ys) #:y-label "Number of cells"
              (list
               (apply render6 (with-x-coords 80 ys))
               (map (lambda (x)
                      (vrule x 0 29 #:width 1 #:style 'long-dash))
                    (range -0.5 70 10))
               (map (lambda (x l)
                      (point-label (vector x 29) l #:size 13 #:point-size 0))
                    (range 0 80 10)
                    (let ((s "Sequence") (o "   Object"))
                      (list s o s o s s o s)))))))]
        [("H3b" "H3c")
         (parameterize
             ([plot-width  200]
              [plot-height 500]
              [plot-x-ticks (linear-ticks #:number 3 #:divisors '(2))]
              [plot-y-ticks (linear-ticks #:number 9 #:divisors '(2))]
              [plot-tick-size 0]
              [plot-x-far-axis? #f]
              [plot-y-far-axis? #f])
           (let ((H3b (string=? figure "H3b"))
                 (H3c (string=? figure "H3c")))
             (display (plot
                       #:title "   Column 1"
                       #:x-min 0 #:x-max 15.3 #:x-label (if H3b "Number of sensations" "")
                       #:y-min 0 #:y-max 4110 #:y-label "Neuron #"
                       (renderH3 L2r L2a L2ac)))
             (when H3c
               (parameterize ([plot-y-axis? #f] [plot-width 150])
                 (display (plot
                           #:title "   Column 2"
                           #:x-min 0 #:x-max 15.3 #:x-label "Number of sensations"
                           #:y-min 0 #:y-max 4096 #:y-label #f
                           (renderH3 L2r1 L2a1 L2a1c)))
                 (display (plot
                           #:title "   Column 3"
                           #:x-min 0 #:x-max 15.3 #:x-label ""
                           #:y-min 0 #:y-max 4096 #:y-label #f
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
              (if H3b (list L2ac) (list L2ac L2a1c L2a2c))
              (if H3b (list L2r) (list L2r L2r1 L2r2)))))]
        ))))

(let loop ()
  (sync (filesystem-change-evt "combined_sequences.data"))
  (sleep 0.1)
  (with-input-from-file "combined_sequences.data"
    (lambda () (plot-cs (read))))
  (loop))

