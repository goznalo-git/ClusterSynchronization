#! /bin/bash

echo > nohup.out

# Enabling multithreading 
export JULIA_NUM_THREADS=40

# Number of simulations
nsims=200

for i in `seq 1 $nsims`; do

    julia ClusterSynchErrors.jl

done
