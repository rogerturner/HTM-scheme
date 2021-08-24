#!chezscheme

;; === HTM-scheme Untangling Sequences Runner  (C) 2019-2021 Roger Turner. ===
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Based on code from Numenta Platform for Intelligent Computing (NuPIC) ;;
  ;; which is Copyright (C) 2017, Numenta, Inc.                            ;;
  ;;                                                                       ;;
  ;; This program is free software: you can redistribute it and/or modify  ;;
  ;; it under the terms of the GNU Affero Public License version 3 as      ;;
  ;; published by the Free Software Foundation.                            ;;
  ;;                                                                       ;;
  ;; This program is distributed in the hope that it will be useful,       ;;
  ;; but WITHOUT ANY WARRANTY; without even the implied warranty of        ;;
  ;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                  ;;
  ;; See the GNU Affero Public License for more details.                   ;;
  ;;                                                                       ;;
  ;; You should have received a copy of the GNU Affero Public License      ;;
  ;; along with this program.  If not, see http://www.gnu.org/licenses.    ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(import
  (chezscheme))

(parameterize (
    (library-directories '(("." . "object-files") ("../../../" . "../../../object-files")))
    (compile-imported-libraries #t)
    (enable-cross-library-optimization #t)
    (optimize-level 3)
    (debug-level 0)
    #;(generate-interrupt-trap #f)
    (generate-inspector-information #f)
    (#%$optimize-closures #t)
    (#%$track-dynamic-closure-counts #f)
    (cp0-effort-limit 20000)             ;; (default is 200)
    (cp0-score-limit   2000)             ;; (default is 20)
    #;(cp0-outer-unroll-limit 1)         ;; (default is 0) [1 slower]
    (generate-wpo-files #t)
    (undefined-variable-warnings #t)
    (compile-file-message #f))
  (compile-program "untangling_sequences.ss" "object-files/untangling_sequences.so")
  (let ([missing (compile-whole-program "object-files/untangling_sequences.wpo" "untangling_sequences.wp")])
    (if (null? missing)
      (let ([command "scheme --program untangling_sequences.wp 6"])
        (display command) (newline)
        (system command)
        (exit))
      (display missing))))
