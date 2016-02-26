# rsa_project
Analysis of gambling fmri data

# file directory

 * `Untitled0.ipynb` inspecting the gambling data
 
 * `a4bres.rds` results of analysis4.R
 * `temp.RData` Temporary file (obv.)  Currently: results for analysis 6, in progress.
 * `temp_pvals.rds` results of analysis 4??
 * `temp_script1.R` Not really a temporary script: used for producing files `roi/data.rds`, `roi/cluster.rds` from Poldrack's 3d data and the clustering results in the `roi/` folder.
 
 * `a6source.R` new bootstrap functions developed for analysis 6 in order to fit subjects individually
 * `analysis1.R` Incomplete.  Attempt to apply ANOVA to voxels/copes.
 * `analysis2.R` Classification of copes.  Finding which ROIs have signal.  Pooling across subjects boosts signal.  DIMENSION REDUCTION.
 * `analysis3.R` Computing distance-feature matrices (without testing, without ordination.)
 * `analysis4_sim_diagnosis.R` At first we got wrong p-values in analysis 4, this simulation study helped lead to realizing that the heterogeneity of covariances was the root of the issue.
 * `analysis4.R` first time we applied bootstrap to get p-values for metric equivalence
 * `analysis5.R` first use of Procrustes
 * `analysis4b.R` update analysis4 by using unbiased test stat
 * `analysis6.R` contrast with analysis 4: pool subjects first then compute test stat
 * `fun_plot.R` cool 3d plot! totally unrelated to the project though...
 * `graphics_a5.R` Plots which try to align MDS distance matrices to the natural parameter distance grid. (Hard to explain..)
 * `indscal_source.R` Find latent parameterizations from distances (when you don't have parametrically generated stimuli...)
 * `instructions.txt` This was an email written by Russell Poldrack on how to get the data from TACC and the experimental parameters.
 * `make_doppelganger_B.R` See `make_doppelganger.R`.  This one tries to imitate the data after it has already been reduced in dimensionality.
 * `make_doppelganger.R` Fits a model to the data in order to generate a synthetic dataset with similar properties.
 * `prepare_dim_reduced.R` An important source file used by many analyses.  Further processes the data matrix, including PCA.
 * `rsa_boot_source.R` A file copied over from github/snarles/fmri then modified.  Contains code for bootstrap-based tests of distance equivalence!
 * `script_output.txt` Output of scripts used to extract ROIS from parametric map data.
 * `source1.R` some functions for working with the nifti 3d brain scans.

## doppel/

Includes synthetic datasets.  (see `make_goppelganger.R`).

## roi/

Various files and scripts involved in extracting ROIs from the parametric brain map.
Also contains the data matrix used in many analyses.

## theory/

Some simulations, mostly about factor analysis (indscal.)