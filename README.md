# TIanalysis
Tuned Inhibition analysis of real fMRI and behavioral data AND simulated fMRI and behavioral data.

TI_LLR_analysis.m is the main file that has 6 components:
 1. Creates the stim (S) parameter that is used in model simulation. It is fitting the stims for the 4 stimulus combinations to the mean d' of all 20 real human subjects. These stims will be used for simulating the LR runs, which is where task behavior is collected.
 
 2. Runs the simulation for LR runs and outputs simulated BOLD, response, and confidence reports. AUC is also calculated for confidence from the x-model, delta-model, and decision (as a potential control)
 2.1 Plots the behavior

 3. Runs the simulation for MI runs and outputs simulated BOLD. The stim for is first manipulated by the same transformation that brought real LR stim combinations to the fitted S for each of the 4 conditions. Inhibition index is also calculated here.

 4. Calculates linear fits (and Pearson, Spearman if needed) of AUC and inhibition index from model simulations (x- and delta- models)

 5. Gathers the correlations from the real data.

 6. Calculates the log likelihood ratio (LLR) of probability of observing the real data given the x-model / the probability of observing the real data given the delta-model.