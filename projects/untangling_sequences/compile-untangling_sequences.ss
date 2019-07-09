(import (except (chezscheme) add1 make-list random reset))

(parameterize (
    (library-directories '(("." . "object-files")))
    (compile-imported-libraries #t)
    (generate-interrupt-trap #f)
    (generate-inspector-information #f)
    (#%$optimize-closures #t)
    (#%$track-dynamic-closure-counts #f)
    (optimize-level 3)
    (cp0-effort-limit 10000)
    (cp0-score-limit 1000)
    (generate-wpo-files #t))
  (compile-program "HTM-scheme/projects/untangling_sequences/untangling_sequences.ss" "object-files/untangling_sequences.so")
  (compile-whole-program "object-files/untangling_sequences.wpo" "untangling_sequences.wp"))
  
(exit)
