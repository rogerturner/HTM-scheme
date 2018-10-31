### HTM algorithms in Scheme

Scheme translations of some [Numenta](https://numenta.com) [HTM](https://numenta.org) (Hierarchical Temporal Memory) algorithms and experiments.

Translated from Numenta [htmresearch](https://github.com/numenta/htmresearch) and [nupic](https://github.com/numenta/nupic).
Directory structure, file names, and source organization echo the htmresearch repository; organization of the code by function, and function and variable names, parallel Numenta code where possible: scheme and Numenta python code can be read side-by-side. Algorithms and project computation code are standard R6RS and run with [Chez Scheme](https://github.com/cisco/ChezScheme) without external libraries. Plotting uses [Racket](http://racket-lang.org) packages; example output from the combined_sequences experiment (replicating figure 6 in [Numenta paper](http://dx.doi.org/10.1101/190678)):


![Figure 6](https://raw.githubusercontent.com/rogerturner/HTM-scheme/master/projects/combined_sequences/Figure%206.png)


 Current algorithms in HTM-scheme/HTM-scheme/algorithms are R6RS libraries named with the path from HTM-scheme root, so can be imported to a top level program in the directory enclosing HTM-scheme.
algorithms/htm_concept.ss includes brief notes on the main data structures used, which are aimed at minimising space requirements.

Archived standalone-spatial-pooler.ss and standalone-temporal-memory.ss are older self-contained Scheme top-level programs including (*obsolete*) library code, some tests, and hello_sp, sp_tutorial, and hello_tm examples. They do not depend on other libraries, so to try these just install [Racket](http://racket-lang.org) or [Chez Scheme](https://github.com/cisco/ChezScheme), Open file in [DrRacket](https://docs.racket-lang.org/drracket/interface-essentials.html), Run, then enter (hello-sp), (sp-tutorial), or (hello-tm), or in Chez Scheme eg:

    $ cd HTM-scheme/archive
    $ scheme
    Chez Scheme Version 9.5.1
    Copyright 1984-2017 Cisco Systems, Inc.

    > (load "standalone-spatial-pooler.ss")
    > (hello-sp)
    See nupic/examples/sp/hello_sp.py
    Random 1:   0.0% (12 99 145 167 198 212 364 412 516 564 569 640 748 758 773 ... 4030)
    Random 2:  13.4% (50 145 153 214 276 379 399 432 538 574 578 663 673 748 824 ... 4034)
    Random 3:   8.5% (12 42 83 124 297 424 451 613 634 666 784 847 859 860 895 ... 4030)
    Repeat 3: 100.0% (12 42 83 124 297 424 451 613 634 666 784 847 859 860 895 ... 4030)
    Noise .1:  97.6% (12 42 83 124 297 424 451 613 634 666 784 847 859 860 895 ... 4030)
    Noise .3:  87.8% (12 42 124 297 424 451 613 634 666 784 847 859 860 895 918 ... 4030)
    ok

[HTM Theory Reading List](https://github.com/rogerturner/HTM-scheme/wiki/HTM-Theory-Reading-List)
