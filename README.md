### HTM algorithms in Scheme

Scheme translations and extensions of some [Numenta](https://numenta.com) [HTM](https://numenta.com/machine-intelligence-technology/htm) (Hierarchical Temporal Memory) algorithms and neuroscience models.

Core algorithms translated from Numenta [htmresearch](https://github.com/numenta/htmresearch) and [nupic](https://github.com/numenta/nupic) repositories. (These Numenta repositories are in maintenance mode: see the [htm.core community fork](https://github.com/htm-community/htm.core) for some [updates](https://github.com/htm-community/htm.core/tree/master/py/htm/advanced/algorithms))

HTM algorithms and models run in [Chez Scheme](https://github.com/cisco/ChezScheme) (standard R6RS Scheme with minor extensions) and are written in "plain" Scheme (ie using only the Scheme forms of common programming language constructs from mainstream PLs). No external libraries or other dependencies are used.

Figure rendering (see examples below) runs in [Racket](http://racket-lang.org) using its [Graph Plotting](https://docs.racket-lang.org/plot/index.html) library.

Example output from the combined_sequences experiment (replicating figure 6 in the Numenta paper [Untangling Sequences- Behavior vs. External Causes](http://dx.doi.org/10.1101/190678)):


![Figure 6](https://raw.githubusercontent.com/rogerturner/HTM-scheme/master/projects/combined_sequences/Figure%206.png)


The HTM-scheme combined_sequences project can also reproduce some of the experiments from [A Theory of How Columns in the Neocortex Enable Learning the Structure of the World](http://dx.doi.org/10.3389/fncir.2017.00081):

#### Fig3(B) One cortical column &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; \(C) Three cortical columns


![Figure 3B/C](https://raw.githubusercontent.com/rogerturner/HTM-scheme/master/projects/combined_sequences/Figure%20H3b+c.png)



### Source organization

Directory structure, file names, and source file structure echo the htmresearch layout; organization of the code by function, and function and variable names, parallel Numenta code where possible: Scheme and Numenta python code can be read side-by-side.

[frameworks/htm-concept.ss](https://github.com/rogerturner/HTM-scheme/blob/master/frameworks/htm-concept.ss) includes brief notes on the main data structures used, which are aimed at minimising space requirements.

Algorithms are implemented as Scheme libraries, named as (directory-name library-name), so
the Chez Scheme (library-directories ...) form should include the path to HTM-scheme.

(Current algorithm implementations may be incompatible with older projects: see CRANKY.md)

Archived standalone-spatial-pooler.ss and standalone-temporal-memory.ss are older self-contained Scheme top-level programs including _(**obsolete**)_ algorithm code, some tests, and hello_sp, sp_tutorial, and hello_tm examples. They do not depend on other libraries, so to try these just install [Racket](http://racket-lang.org) or [Chez Scheme](https://github.com/cisco/ChezScheme), Open file in [DrRacket](https://docs.racket-lang.org/drracket/interface-essentials.html), Run, then enter (hello-sp), (sp-tutorial), or (hello-tm) -- or in Chez Scheme:

    $ cd HTM-scheme/archive
    $ scheme
    Chez Scheme Version 9.5.2
    Copyright 1984-2019 Cisco Systems, Inc.

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

For an introduction to development in Scheme, the EdX course "How to Code" is recommended:

    "For students who have already programmed, this module may seem too easy – 
    but be careful! Be sure to learn **this** programming language.
    While this language forms the conceptual core of nearly every programming
    language you might have used, it is also different in important ways.
    Take the time to go through this material carefully. Working through the
    videos and practice materials for this module should take approximately
    5-8 hours of dedicated time to complete."
[EdX How to Code Module 1a Overview](https://learning.edx.org/course/course-v1:UBCx+HtC1x+2T2017/) (course consists of 10 modules)

[HTM Theory Reading List](https://github.com/rogerturner/HTM-scheme/wiki/HTM-Theory-Reading-List)
