**Attention***: Handling this repository may require the Git LFS extension to be installed on your system.

This is the git repository containing all code needed to reproduce most computationally generated Figures presented in "Mesenchymal stromal cell-derived septoclasts resorb cartilage during developmental ossification and fracture healing" by Sivaraj et al. published in Nature Communications, 2022. This does not include Figures related to the Fracture Experiments presented in the second half of the Paper.



The repository, like the analysis performed, is split into three parts and their associated folders:



1. **Preprocessing**: This folder Contains two sub-folders, one for each of the general datasets used in the Analysis. These folders contain a snakamake pipeline each to download the appropriate genome and process the read files placed into */data/reads/ (they may need to be named approprietly to match the presets of the Pipeline). The Pipeline can then be run by adapting and executing the "run_scRNA.sh" in the respective folder. These Pipelines will output read matrices, which can be read into R in the next step. These workflows come with H2B-GFP sequences downloaded from AddGene at https://www.addgene.org/11680/sequences/. Also included is the official 10X Chromium V3 whitelist. This file was copied from [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger).
2. **Processing**: This Folder contains an RMarkdown Workflow and an accessory script containing helper functions. Simply running the main RMarkdown script (in for example RStudio) from the main project folder creates all files necessary for the creation of all main Figures in the following step.
3. **Figure_Creation**: This folder contains another RMarkdown Workflow. Running this script from the Projects main Folder will create a sub-folder "Figures" in "Figure_Creation", containing all main and supplementary figures.

Be aware, that some methods (such as UMAP) used in this Workflow have components that generate stochastic output. This was not properly accounted for when initially preparing the Figures. Therefore, the Figures produced by running this script locally might differ slightly from the ones presented in the publication. To exactly reproduce these Figures, the original R objects from which the Figures in the publication were created, will be shared upon request.





Estimated Runtime of entire Workflow, using max. 40 cores: 3h 



