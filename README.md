# Cluster Synchronization Errors in Julia

This is the repository hosting some basic scripts and data used in the papers
- "Endowing networks with desired symmetries and modular behavior". P. Khanra, S. Ghosh, D. Aleja, K. Alfaro-Bittner, G. Contreras-Aso, R. Criado, M. Romance, S. Boccaletti, P. Pal, and C. Hens; [https://arxiv.org/abs/2302.10548](https://arxiv.org/abs/2302.10548).
- "The transition to synchronization of networked systems". A. Bayani, F. Nazarimehr, S. Jafari, K. Kovalenko, G. Contreras-Aso, K. Alfaro-Bittner, R.J. Sánchez-Garcı́a, and S. Boccaletti.

## General information

This repository is not an exact replica of the actual project. By this we mean that the scripts can't be run "as is", and the data shown here is not the one directly obtained by running the scripts. Instead, the data files in `Data_all_simulations` contain the already pre-processed (essentially averaged) values of the synchronization errors per cluster of each simulation, namely:
- `smallweighted_lorenz.json`: N=10 fully connected, weighted network, with Lorenz dynamics and coupling on the X coordinates.
- `smallweighted_rossler.json`: N=10 fully connected, weighted network, with Rossler dynamics and coupling on the Y coordinates.
- `synth_network_1000.json`: N=1.000 synthetic network, with Rossler dynamics and coupling on the Y coordinates.
- `synth_network_10000.json`: N=10.000 synthetic network, with Rossler dynamics and coupling on the Y coordinates.
- `real_network_powergrid.json`: USA Power grid network, with Rossler dynamics and coupling on the Y coordinates. a,b,c fixed.
- `powergrid_variable-a01.json`: USA Power grid network, with Rossler dynamics and coupling on the Y coordinates. 10% heterogeneity in a.
- `powergrid_variable-b01.json`: USA Power grid network, with Rossler dynamics and coupling on the Y coordinates. 10% heterogeneity in b.
- `powergrid_variable-c01.json`: USA Power grid network, with Rossler dynamics and coupling on the Y coordinates. 10% heterogeneity in c.

As for the scripts `ClusterSynchErrors-?????.jl`, both are written in Julia (keep in mind it was our first encounter with such programming language, and learn the bare minimum to carry out this project, so it could be qualified as "spaguetti code"), using some basic scientific libraries (DifferentialEquations, MAT, LinearAlgebra, Statistics...) and making use of multithreading, parallelizing the solver with a different value of the coupling per thread. 

To automatize running the Julia script over and over again, the `batch_senders.sh` bash file can be run using
```bash
nohup ./batch_senders.sh &
```