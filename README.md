# Cluster Synchronization Errors in Julia

This is the repository hosting some basic scripts and data used in the papers
- "Endowing networks with desired symmetries and modular behavior". P. Khanra, S. Ghosh, D. Aleja, K. Alfaro-Bittner, G. Contreras-Aso, R. Criado, M. Romance, S. Boccaletti, P. Pal, and C. Hens; [Phys. Rev. E 108, 054309 (2023)](https://doi.org/10.1103/PhysRevE.108.054309).
- "The transition to synchronization of networked systems". A. Bayani, F. Nazarimehr, S. Jafari, K. Kovalenko, G. Contreras-Aso, K. Alfaro-Bittner, R.J. Sánchez-Garcı́a, and S. Boccaletti; [Nat. Comm. 15, 4955 (2024)]([https://doi.org/10.48550/arXiv.2303.08668](https://doi.org/10.1038/s41467-024-48203-6))

## General information

The repository is structured as follows:

- Folder `SynchronizationErrors` contains the scripts necessary for computing the cluster synchronization errors of a given dynamical system over a given network. They are written in Julia using some basic scientific libraries (DifferentialEquations, MAT, LinearAlgebra, Statistics...) and making use of multithreading, parallelizing the solver with a different value of the coupling per thread.
  - The `ClusterSynchErrors-*.jl` scripts can be modified with the desired dynamical system and network, and then run.
  - The `class*dynamics.jl` contains a collection of possible dynamical systems to consider, already formatted appropriately.
  - To automatize running the Julia script over and over again, the `batch_senders.sh` bash file can be run using
  ```bash
  nohup ./batch_senders.sh &
  ```

- Folder `Input` contains adjacency and/or Laplacian matrices of some networks, as well as initial conditions (close to the respective attractors).

- The data files in `Data_all_simulations` contain the already pre-processed (averaged) values of the synchronization errors per cluster of each simulation shown in the aforementioned papers.

- Folder `LyapunovExpCalculation` contains scripts which can be run to compute the Master Stability Function (maximum Lyapunov exponent) of a given dynamical system, along with examples of them. 
