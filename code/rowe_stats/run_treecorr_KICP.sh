#!/bin/bash                                                                                                                                                   
echo "RUNNING"
echo "slurm procid = " $SLURM_PROCID
echo "slurm ntasks = " $SLURM_NTASKS
echo "slurm cpus per task = " $SLURM_CPUS_PER_TASK

#Compute rho stats measurement                                                                                                                           
cmd="python ./run_treecorr.py"
date
echo $cmd
$cmd
date
