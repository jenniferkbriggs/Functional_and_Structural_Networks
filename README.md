# Functional_and_Structural_Networks
This is all the code used in the paper: *Briggs, Jennifer K., et al. "Beta-cell intrinsic dynamics rather than gap junction structure dictates subpopulations in the islet functional network." Elife 12 (2023): e83147.*

This directory is used for running network analysis on time courses. Particularly, this code was developed for calcium and electrical signals in the Islets of Langerhans but may be used for other systems.

Managed by Jennifer Briggs 2022.


## Getting started
In Briggs et al., Elife 2023, we study the relationships between structural and functional networks in the islet of Langerhans. We use two simulations of the islet and experimental data. The Cha-Noma fast oscillation model of the islet can be found in folder: *FastOscillation_sim* (this code is not fully commented). 

The easiest way to explore analysis for all *simulated data* is to use the *RunAllAnalysis.m* script. This script will go through analysis for Figures 1, 4, 5, and 7 in the paper. 

## Example Data
This directory has example data from the manuscript that can be used to learn about the code. Due to space limitations on github, we uploaded two seeds of the small 300 cell fast oscillation model and one seed of the downsampled coupled slow oscillation model (IOM). These are not the exact same data used in the paper. For this reason, the results may differ. To simulate the 1000 cell fast oscillation model, use code in **FastOscillation_Sim**. To simulate the coupled slow oscillation model (IOM), contact Dr Isabella Marinelli (i.marinelli@bham.ac.uk). We also include one csv of example experimental data. 


