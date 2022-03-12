### Simulation of MR results to choose mediation method

The repo contains simulation code to help decide what mediation method to use for a given set of GWAS/MR data/results. This was a supplementary analysis in [this](https://github.com/mvab/mendelian-randomization-breast-cancer) project.


`submit_script.sh` runs `process.R` (uses functions from `original_code.R`). Then, the simulation output is summarised in `review_output.R`.

To submit simulation script on BC3: `qsub submit_script.sh` 

 - set walltime if > 500 iters `qsub -l walltime=05:00:00`
 - specify cores, iters, mediator to process in .sh
 - to set mediator row: `qsub submit_script.sh -F "2"` 



[![DOI](https://zenodo.org/badge/314559039.svg)](https://zenodo.org/badge/latestdoi/314559039)

