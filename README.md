### HTM Spatial Pooler and Temporal Memory algorithms in Scheme

The current Scheme versions of these algorithms are the libraries lib-sp.ss and lib-tm.ss,
and are used by programs in the examples folder (configured for Chez Scheme).

Translated from spatial_pooler.py, temporal_memory.py, hello_sp.py, sp_tutorial.py, hello_tm.py, and tm-high-order.py in 
[Numenta NuPIC](https://github.com/numenta/nupic): see comments there for details.

spatial-pooler.ss and temporal-memory.ss are self-contained Scheme R<sup>6</sup>RS top-level programs including (older) library code, some tests, and hello_sp, sp_tutorial, and hello_tm examples.

To get started with these just install [Racket](http://racket-lang.org) or [Chez Scheme](https://github.com/cisco/ChezScheme),

Open file in [DrRacket](https://docs.racket-lang.org/drracket/interface-essentials.html) 
and Run, or (load "file.ss") in Chez Scheme repl, then (hello-sp), (sp-tutorial), or (hello-tm)
