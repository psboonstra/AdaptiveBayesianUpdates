# Adaptive Bayesian Updates

This is a companion repository for the `R` package `adaptBayes` available here: https://github.com/umich-biostatistics/adaptBayes

After installing and loading that package using the following commands,

```r
library(devtools)
# may take some time:
install_github('umich-biostatistics/adaptBayes') 

library(adaptBayes)
```

then use the scripts in this repository to reproduce 
the simulations studies in Boonstra and Barbaro

## Further details

In more detail, there are six files included in this repository (in addition to this README): one text file (ending in <samp>.txt</samp>) and five <samp>R</samp> scripts (ending in  <samp>.R</samp>). The simulation studies reported in Boonstra and Barbaro were run using commit 22

### Text file
<samp>run_abu_sims.txt</samp> is the script for submitting parallel runs of <samp>run_aub_sims.R</samp> (described below) to a cluster that is running SLURM. The following command run in the terminal will do this:

<code> sbatch run_abu_sims.txt </code>

The script assumes that you are working in your home directory. Edit the script as necessary  

### <samp>R</samp> files

<samp>exemplar.R</samp> creates a single simulated dataset and walks through analyzing these data using the various adaptive priors

<samp>run_abu_sims.R</samp> is the script to conduct the large-scale simulation study described in the manuscript. On a local machine, the user may choose a specific <samp>array_id</samp> (as described in this script's documentation) and run the code locally on his/her machine. On a cluster running SLURM, the user can use this script to submit multiple jobs simultaneously (as described above). 

<samp>functions_simulation.R</samp> provides the simulation functions. NOTE: this only contains the extra functions needed to simulate the data; the methods are contained in the `adaptBayes` package

<samp>generate_params.R</samp> constructs inputs for running the simulation study. As described in the script's documentation and the language below, these inputs can be overwritten by the user

<samp>make_figures.R</samp> gives the code to create the figures and tables in the manuscript and supplementary material reporting on the simulation study. 

### Current Suggested Citation

Boonstra, Philip S. and Barbaro, Ryan P., "Incorporating Historical Models
with Adaptive Bayesian Updates" (2018) *Biostatistics* 
https://doi.org/10.1093/biostatistics/kxy053

<a href="https://biostats.bepress.com/umichbiostat/paper124">Authors' Copy </a>

DOI for this repository:

[![DOI](https://zenodo.org/badge/140338593.svg)](https://zenodo.org/badge/latestdoi/140338593)
