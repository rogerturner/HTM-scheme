(import (except (chezscheme) add1 make-list random reset))

(parameterize (
    (library-directories '(("." . "object-files")))
    (compile-imported-libraries #t)
    (generate-interrupt-trap #f)
    (generate-inspector-information #f)
    (#%$optimize-closures #t)
    (#%$track-dynamic-closure-counts #f)
    (optimize-level 3)
    (cp0-effort-limit 10000)             ;; (default is 200)
    (cp0-score-limit 1000)               ;; (default is 20)
    (enable-cross-library-optimization #t)
    (generate-wpo-files #t))
  (compile-program "HTM-scheme/projects/KropffTreves2008_reproduction/grid_cell_demo2.ss" "object-files/grid_cell_demo2.so")
  (compile-whole-program "object-files/grid_cell_demo2.wpo" "grid_cell_demo2.wp"))
  
(exit)
