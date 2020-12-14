Simulation of MR results to choose mediation method


To submit simulation script on BC3: `qsub submit_script.sh` 

 - set walltime if > 500 iters `qsub -l walltime=05:00:00`
 - specify cores, iters, mediator to process in .sh
 - to set mediator row: `qsub submit_script.sh -F "2"` 

`submit_script.sh` runs `process.R` (uses functions from `original_code.R`). Then, the simulation output is in  summaries in `review_output.R`.