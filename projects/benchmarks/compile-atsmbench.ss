(import (chezscheme))

(parameterize (
    (library-directories '(("." . "object-files")))
    (generate-interrupt-trap #f)
    (compile-imported-libraries #t)
    (optimize-level 3)
    (cp0-effort-limit 10000)
    (cp0-score-limit 1000)
    (cp0-outer-unroll-limit 50)
    (generate-wpo-files #t))
  (compile-program "HTM-scheme/projects/benchmarks/atsmbench.ss")
  (compile-whole-program "HTM-scheme/projects/benchmarks/atsmbench.wpo" "atsmbench.wp"))
  
(exit)
