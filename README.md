### HTM algorithms in Scheme

[Chez Scheme](https://cisco.github.io/ChezScheme/) translations of some [Numenta](https://numenta.com) [HTM](https://numenta.org) (Hierarchical Temporal Memory) algorithms.

Translated from spatial_pooler.py, temporal_memory.py, hello_sp.py, sp_tutorial.py, hello_tm.py, and tm-high-order.py in [Numenta NuPIC](https://github.com/numenta/nupic) and ApicalTiebreakTemporalMemory.cpp in [htmresearch-core](https://github.com/numenta/htmresearch-core): see comments there for details.


standalone-spatial-pooler.ss and standalone-temporal-memory.ss are older self-contained Scheme R<sup>6</sup>RS top-level programs including (*obsolete*) library code, some tests, and hello_sp, sp_tutorial, and hello_tm examples. To get started with these install [Racket](http://racket-lang.org) or [Chez Scheme](https://github.com/cisco/ChezScheme), Open file in [DrRacket](https://docs.racket-lang.org/drracket/interface-essentials.html) and Run, or (load "file.ss") in Chez Scheme repl, then (hello-sp), (sp-tutorial), or (hello-tm)
