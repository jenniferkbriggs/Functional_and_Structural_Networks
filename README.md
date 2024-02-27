# Functional_and_Structural_Networks
This is all the code used in the paper: *Briggs, Jennifer K., et al. "Beta-cell intrinsic dynamics rather than gap junction structure dictates subpopulations in the islet functional network." Elife 12 (2023): e83147.*

This directory is used for running network analysis on time courses. Particularly, this code was developed for calcium and electrical signals in the Islets of Langerhans but may be used for other systems.

Managed by Jennifer Briggs 2022. If you have particular questions, feel free to reach out: Jennifer.briggs@cuanschutz.edu



## Getting started
In Briggs et al., Elife 2023, we study the relationships between structural and functional networks in the islet of Langerhans. We use two simulations of the islet and experimental data. The Cha-Noma fast oscillation model of the islet can be found in folder: **FastOscillation_sim** (this code is not fully commented). 

The easiest way to explore analysis for all **simulated data** is to use the **RunAnalysis_simulation.m** script in the Analyze Data folder. This script will go through analysis for Figures 1, 4, 5, and 7 in the paper. 

## Example Data
This directory has example data from the manuscript that can be used to learn about the code. Due to space limitations on github, we uploaded two seeds of the small 300 cell fast oscillation model and one seed of the downsampled coupled slow oscillation model (IOM). These are not the exact same data used in the paper. For this reason, the results may differ. To simulate the 1000 cell fast oscillation model, use code in **FastOscillation_Sim**. To simulate the coupled slow oscillation model (IOM), contact Dr Isabella Marinelli (i.marinelli@bham.ac.uk). We also include one csv of example experimental data. (You must unzip this directory before using it) 


## Analyze Data
This directory contains all of the scripts used to analyze data. 

**RunAnalysis_simulation.m** will run simulation analyses for the fast simulation corresponding to Figures 1, 4, 5, and 7 in the paper. 

**AnalyzeExperimental.m** was used to create Figure 3 by analyzing experimental calcium traces and FRAP data. Note that the calcium traces are assumed to be in .csv format (see extract data). FRAP is also in .csv. with rows containing the rate of recovery of the corresponding cell. We refer the reader to Farnsworth, Nikki L., et al. The Journal of physiology 592.20 (2014): 4431-4446. for how to analyze FRAP data. 

**AnalyzeIOM.m** is used to analyze the IOM data.

## Extract Data
This repository contains code to convert imaging files into calcium csv (such as in ExampleData/experimental/exampleCA.csv). 

When conducting network analysis, it is very important that the pixels which we consider part of cell (i) are actually part of cell (i) and not cell (j), for example. **STanalysis.m** helps us ensure that the pixels within a single cell mask are truly part of that cell.

**ExtractImages.m** is a script used to convert .lsm/.czi files into calcium traces. This code calls **STanalysis.mat** and can be used to understand the CellMask variable loaded into STanalysis.m. However, because most people use ImageJ or some other software to extract the cell files, this code is not commented in detail and is not currently compatible with all computers and file structures.  

## Visualization
These are scripts that I made to visualize the network. They are not currently well commented. 
