# R code for the simulation study reported in Boonstra and Barbaro (2018)

This is a companion repository for the `R` package `adaptBayes` available at https://github.com/umich-biostatistics/adaptBayes. The `adaptBayes` package
contains the code for actually fitting an adaptive prior in a Bayesian GLM. 
*This* repository only contains the additional code needed for conducting the 
simulation study in Boonstra and Barbaro (2018). Before you can use this 
repository, you need to first install and load the `adaptBayes` package using 
the following commands:

```r
if(!require(adaptBayes)) {
  library(devtools)
  # if installation is necessary, compiling everything will take a few minutes
  install_github('umich-biostatistics/adaptBayes') 
} 
```

After installing the package, you can use the scripts in this repository to 
reproduce the simulations studies in the manuscript. See also `vignette.pdf` 
for a vignette describing the typical usage of the adaptive Bayesian priors on
a simulated dataset

### Further details

In more detail, there are six files included in this repository (in addition to 
this README and `vignette.pdf`): one text file (<samp>run_abu_sims.txt</samp>) 
and five <samp>R</samp> scripts (ending in  <samp>.R</samp>). The simulation
studies reported in Boonstra and Barbaro were run using commit 22

#### Text file
<samp>run_abu_sims.txt</samp> is the script for submitting parallel runs of
<samp>run_aub_sims.R</samp> (described below) to a cluster that is running
SLURM. The following command run in terminal will do this:

<code> sbatch run_abu_sims.txt </code>

The script assumes that you want all of the results to be put in your home 
directory (which you probably don't). Edit the script as needed  

#### <samp>R</samp> files

<samp>vignette.R</samp> creates a single simulated dataset and walks through 
analyzing these data using the various adaptive priors. It can also be `knit` 
to create a copy of `vignette.pdf` (it will take a few minutes to knit)

<samp>run_abu_sims.R</samp> is the script to conduct the large-scale simulation
study described in the manuscript. On a local machine, the user may choose a
specific <samp>array_id</samp> (as described in this script's documentation) 
and run the code locally on his/her machine. On a cluster running SLURM, the
user can use this script to submit multiple jobs simultaneously (as described 
in the description of <samp>run_abu_sims.txt</samp> above). 

<samp>functions_simulation.R</samp> provides the simulation functions. NOTE: 
this only contains the extra functions needed to simulate the data; the methods
are contained in the `adaptBayes` package

<samp>generate_params.R</samp> constructs inputs for running the simulation
study. As described in the script's documentation and the language below, these
inputs can be overwritten by the user

<samp>make_figures.R</samp> gives the code to create the figures and tables in 
the manuscript and supplementary material reporting on the simulation study. 

#### Current Suggested Citation

Boonstra, Philip S. and Barbaro, Ryan P., "Incorporating Historical Models with
Adaptive Bayesian Updates" (2018) *Biostatistics* 
https://doi.org/10.1093/biostatistics/kxy053

<a href="https://biostats.bepress.com/umichbiostat/paper124">Authors' Copy </a>

DOI for this repository:

[![DOI](https://zenodo.org/badge/140338593.svg)](https://zenodo.org/badge/latestdoi/140338593)
