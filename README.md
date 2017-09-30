### HTM Spatial Pooler and Temporal Memory algorithms in Scheme

The current Scheme versions of these algorithms are the libraries htm-sp.ss, htm-tm.ss, and htm-attm.ss, and are used by programs in the examples folder (configured for Chez Scheme).

Translated from spatial_pooler.py, temporal_memory.py, ApicalTiebreakTemporalMemory.cpp, hello_sp.py, sp_tutorial.py, hello_tm.py, and tm-high-order.py in [Numenta NuPIC](https://github.com/numenta/nupic): see comments there for details.


spatial-pooler.ss and temporal-memory.ss are older self-contained Scheme R<sup>6</sup>RS top-level programs including (*obsolete*) library code, some tests, and hello_sp, sp_tutorial, and hello_tm examples. To get started with these install [Racket](http://racket-lang.org) or [Chez Scheme](https://github.com/cisco/ChezScheme), Open file in [DrRacket](https://docs.racket-lang.org/drracket/interface-essentials.html) and Run, or (load "file.ss") in Chez Scheme repl, then (hello-sp), (sp-tutorial), or (hello-tm)
