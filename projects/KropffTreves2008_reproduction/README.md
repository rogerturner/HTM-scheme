# Reproduction<sup>2</sup> of Kropff & Treves, 2008
### Grid cells using HTM-scheme  

This is a reproduction study of:

	The emergence of grid cells: intelligent design or just adaptation? 
	Emilio Kropff and Alessandro Treves, 2008
	
translated from the [NUPIC based reproduction](https://github.com/ctrl-z-9000-times/KropffTreves2008_reproduction) by David McDougall into HTM-scheme.

#### This README is a lightly edited abridgement of https://github.com/ctrl-z-9000-times/KropffTreves2008_reproduction/blob/master/README.md by David McDougall

## Introduction
Kropff & Treves [1] describe a method of producing grid cells.  It appears to be supported by biology [2].  The gist of it is to make a spatial pooler [3] with two new mechanisms: stability and fatigue.  The spatial pooler's overlap is stabilized by putting it through a low pass filter, which smooths over input fluctuations and removes spurious input fluctuations.  Each cell also has a fatigue which slowly follows its activity, catches up to it and turns it off.  The fatigue is equal to the overlap passed through a low pass filter.  
* The effect of the stability mechanism is to cause the cells (or SP mini-columns) to learn large, contiguous areas of the input.  As the sensory organ moves around the world, the stability mechanism causes the cells to represent large & contiguous areas of the world by forcing cells to react slower than their sensory input is changing.  
* The fatigue mechanism shapes the contiguous areas of the grid cell receptive fields into spheres which are then packed into the environment.  As the sensory organ passes through areas of the world, grid cells get tired and fall behind in the competition for a short while.  

## Methods
Grid cells are implemented as an extension to the Spatial Pooler (SP) algorithm.  The SP is modified in three ways:
1)	Stability and fatigue dynamics are applied to the overlap.
2)	The overlap is divided by the number of connected synapses to each grid cell.
3)	Synapses only learn when either the presynaptic or postsynaptic side changes its activity state.  This filters out duplicate updates on sequential time steps.  This causes the grid cells to only learn when the agent is moving around, a stationary agent's grid cells will not learn.

#### Usage

    $ scheme
    > (load "HTM-scheme/projects/KropffTreves2008_reproduction/compile-grid_cell_demo.ss")
    $ scheme --program grid_cell_demo.wp -tt 200000

    and then run HTM-scheme/projects/generate_plots.rkt in Racket

## Results

Compare these results to Kropff & Treves, 2008.

The model is trained in a 256x256 arena by randomly walking.  The simulated agent moves at a constant speed of 1.4 units per step.  When it reaches one of the arena's boundaries it is turned to face back into the arena.  The model is trained for 200,000 steps, which took 20 minutes.

![Agent's path through arena, 10k steps](Path_10k.png?raw=true "Path, 10k steps")

**Figure 1:** Example of random walk.  This figure contains only 10,000 steps.

---

Each location in the arena is represented by sparse input from a total of 2,500 place cells. 75 to 150 cells are active for each location, for sparsity of 3% - 6%.  The coordinate encoder models these input cells.

![Place Cell Receptive Fields](Input_Receptive_Fields.png?raw=true "Place Cell Receptive Fields")
**Figure 2:** Example place cell receptive fields.  Each plot is for a randomly selected place cell.  The place cell activates when the agent moves into any of the indicated areas of the arena.

---

The model is tested by examining which locations each grid cell activates at (AKA its receptive field).  Learning is disabled while testing.  Figure 3 shows the results of this test performed on an untrained model.  Figure 4 shows the results of this test performed on a trained model.

![Untrained Grid Cell Receptive Fields](Grid_Cell_Receptive_Fields_untrained.png?raw=true "Untrained Grid Cell Receptive Fields")
**Figure 3:** Untrained grid cell receptive fields.  Each box is a randomly selected grid cell.
Notice that the grid cells do not respond to large contiguous areas of the arena.  There are many isolated (non-contiguous) activations.  The small contiguous area are randomly shaped and have fuzzy, ill-defined borders.  

![Trained Grid Cell Receptive Fields](Grid_Cell_Receptive_Fields_trained.png?raw=true "Trained Grid Cell Receptive Fields")
**Figure 4:** Receptive fields of grid cells after training (top row: same cells as figure 3, other rows: first 10 columns).  Notice that the grid cells respond to large contiguous areas.  Many of the receptive fields are approximately round, the same size, and have sharp, well defined borders.  

## References
[1]	The emergence of grid cells: intelligent design or just adaptation? Emilio Kropff and Alessandro Treves, 2018.  DOI 10.1002/hipo.20520

[2]	Tuning of Synaptic Integration in the Medial Entorhinal Cortex to the Organization of Grid Cell Firing Fields, Garden, Dodson, O’Donnell, White, Nolan, 2008.  DOI 10.1016/j.neuron.2008.10.044

[3]	Cui Y, Ahmad S and Hawkins J (2017), The HTM Spatial Pooler—A Neocortical Algorithm for Online Sparse Distributed Coding. Front. Comput. Neurosci. 11:111. doi: 10.3389/fncom.2017.00111
