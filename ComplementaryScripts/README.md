# ecYeast cell-factory

# Instructions

* All scripts are refered relative to their location inside of the subfolder `cellFactory-yeast-GEM/ComplementaryScrips`.

1. Run MATLAB script `eval_fluxDist.m` to evaluate biomass yield, product yield and product rate for each of the +100 chemicals using ecYeastGEM (wild-type). This script also produces a matrix of stoichiometric distances between optimal distributions for maximum production of each chemical.

2. Run MATLAB script `method/find_gene_targets.m` to run the `robust_ecFSEOF` method on all the models available in this repository. This will generate lists of gene target predictions for enhanced production of +100 chemicals in *S. cerevisiae*.

3. Run MATLAB scripts `summarize_results.m` and `getTargetsMatrix.m` to obtain a summary of the number of predicted targets for OE, KD and KO for each product, and a representation of this in a discrete matrix form, respectively.

4. Run R script `analyse_results.R`. This script will generate:
 	- a pie plot showing the products classification by chemical family. 
 	- Preliminary cluster analysis of gene targets for all products (hclust heatmaps) and dimensionality reduction of gene targets data (PCA and t-SNE plots)
 	- A heatmap showing euclidian distances between all vectors of targets (for each product). In order to so, gene OEs are assigned a value=2, value=0 for KOs, value = 0.5 for KDs and 1 for unaffected genes.
 	- Boxplots, showing distribution of number of targets, per type, for all products and levels of the multi-step robust ecFSEOF (1. FSEOF, 2. mechanisctically validated individual targets, 3. optimal combination of targets).

 

