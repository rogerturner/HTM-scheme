
*(Adapted from numenta/htmresearch/projects/combined_sequences/README.md)*
## Untangling Sequences: Behavior vs. External Causes

This directory contains Scheme code that runs sensorimotor sequences
in combination with pure temporal sequences, using the HTM-scheme
implementation of Numenta NuPIC algorithms.

The goal is to test whether a single neural mechanism can
automatically discover what parts of a changing sensory stream are due
to movement and which parts are due to external causes, and to learn
predictive models of both types of causes simultaneously using simple
learning rules.

### Abstract (From Numenta [draft paper](http://dx.doi.org/10.1101/190678))

There are two fundamental reasons why sensory inputs to the brain change
over time. Sensory inputs can change due to external factors or they can
change due to our own behavior. Interpreting behavior-generated changes
requires knowledge of how the body is moving, whereas interpreting
externally-generated changes relies solely on the temporal sequence of
input patterns. The sensory signals entering the neocortex change due to
a mixture of both behavior and external factors. The neocortex must
disentangle them but the mechanisms are unknown. In this paper, we show
that a single neural mechanism can learn and recognize both types of
sequences. In the model, cells are driven by feedforward sensory input
and are modulated by contextual input. If the contextual input includes
information derived from efference motor copies, the cells learn
sensorimotor sequences. If the contextual input consists of nearby
cellular activity, the cells learn temporal sequences. Through
simulation we show that a network containing both types of contextual
input automatically separates and learns both types of input patterns.
We review experimental data that suggests the upper layers of cortical
regions contain the anatomical structure required to support this
mechanism.

### Usage
To run one of the experiments using Chez Scheme (threaded 32-bit) and Racket, say the one for Figure 4A:

    $ cd <directory containing HTM-scheme directory>
    $ scheme
    > (load "HTM-scheme/projects/untangling_sequences/untangling_sequences.ss")
    > (run exp4a AH2017)
    
    Then in DrRacket File>Open… projects/generate_plots.rkt and Racket>Run

Note: the results may not be identical to the charts in the paper, due
to changes in the random number generator
and perhaps also algorithm changes.  They should be similar though, and
conclusions and takeaways should be the same.

