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
  (compile-program "HTM-scheme/projects/combined_sequences/combined_sequences.ss")
  (compile-whole-program "HTM-scheme/projects/combined_sequences/combined_sequences.wpo" "combined_sequences.wp"))
  
(exit)
