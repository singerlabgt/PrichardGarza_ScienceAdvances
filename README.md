# PrichardGarza_ScienceAdvances

## Overview
This code was used to generage figures for Prichard, Garza, et. al., published in Science Advances in 2023. This code generates data structures, graphs, and statistical analyses based on measures from mice exposed to flickering auditory or visual stimuli. Tissue analyses including changes in microglia morphology were measured using the Imaris imaging software to reconstruct cells (Figure 1 & Figure 6) and pNFkB colocalization within cells (Figure 5) were examined in visual cotex tissue samples from one hemisphere. The other hemisphere visual cortex was examined for cytokine levels (Figure 3 & Figure 7) as well as RNA single cell data (Figure 2 & Figure 4). Final versions for all figures for the publication were edited in Adobe Illustrator.

## Requirements
Matlab 2020B

.csv files including cytokine data

.csv files including microglia morphology data from Imaris

Violin Plot https://github.com/bastibe/Violinplot-Matlab

## Figure/ Code Overview

### Figures 1 and Figure 6
These figures were generated by using matlab code located in the glia_morphology subfolder to create a data structure of morphology measures including branch length, branch depth, process length, process volume, etc.
These measures were then compared across groups using ANOVAs then passed to Violinplot functions to create graphs. Reconstructed microglia images are stored on NeuroMorph.

![image](https://github.com/singerlabgt/PrichardGarza_ScienceAdvances/assets/57195922/052c7b07-198e-4fb9-a755-c3a7120c73db)

### Figure 5
The colocalization of pNfkB within neurons or microglia were measured using Imaris software to construct surfaces. Matlab colocalization_code was used to divide the amount of pNFkB within cells by the volume of each cell type, then a ratio was created. Graphs were then created to compare pNFkB volume between cells, then statistically compared using an Anova. 

![image](https://github.com/singerlabgt/PrichardGarza_ScienceAdvances/assets/57195922/daf7324b-10f2-41ff-9337-901684e1b2c8)

### Figure 3 and Figure 7 
These were generated using cytokine data in the attached .csv file.  PLSDA graphs and heatmaps were generated using the PLSDAFunction. Permutations scripts were used to randomly permute the sample names 1k times to create a null distribution to show that the results were not due to chance etc. 

![image](https://github.com/singerlabgt/PrichardGarza_ScienceAdvances/assets/57195922/dc4f04e0-281e-4014-85d1-ca808776d794)

The scripts also generate individual cytokine bar plots for each of the specified cytokines. These graphs as well as the heatmap are recolored using the folder othercolor. The final colors for the graphs and marks over individual bars to denote significance values were edited using Adobe Illustrator.

![image](https://github.com/singerlabgt/PrichardGarza_ScienceAdvances/assets/57195922/5e9d04fb-6e47-4aac-b2cd-4ab5e2885c04)

### Figure 2
This figure of immunity gene sets were generated using the code in the Figure 2 Code and Data folder. Data is included in the .csv file within the folder. 

![image](https://github.com/singerlabgt/PrichardGarza_ScienceAdvances/assets/57195922/260ee468-c1b7-4f43-823f-03f7051257af)







