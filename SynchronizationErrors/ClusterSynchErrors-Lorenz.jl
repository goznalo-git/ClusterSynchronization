using DifferentialEquations
using MAT
using SparseArrays
using LinearAlgebra
using Statistics
using NPZ

include("cluster_error_functions.jl")

####################################
########## LOADING INPUTS ##########
####################################

# Load directly the Laplacian matrix
file1 = matread("Input/L_SmallWeighted.mat");
Lij = file1["L"];

#Sparsify the Laplacian matrix
sparseL = sparse(Lij);


# Load the Lorenz initial condition (within the attractor)
file2 = matread("Input/initc_Lorenz.mat");
initc = file2["initc"];


##########################################
########## ESTABLISH PARAMETERS ##########
##########################################


# Parameters and initial conditions
N = size(sparseL, 1);
x0 = initc[1] .+ 0.01*rand(N);
y0 = initc[2] .+ 0.01*rand(N);
z0 = initc[3] .+ 0.01*rand(N);
u0 = [x0 y0 z0];

len = 400 # Number of couplings (lambda) to calculate over 
lambdas = range(0, 8, length = len)


# Check the number of threads available for the parallelization
nth = Threads.nthreads()
println(nth);


##########################################
########## EQUATIONS AND SOLVER ##########
##########################################

# Establish ODE 
function lorenz!(du, u, lambda, t, sigma=10, rho=28, beta=2)
    coupling = similar(x0);
    mul!(coupling, sparseL, u[:,1]);
    du[:,1] = sigma * (u[:,2] .- u[:,1]) .- lambda .* coupling;
    du[:,2] = u[:,1] .* (rho .- u[:,3]) .- u[:,2]; 
    du[:,3] = u[:,1] .* u[:,2] .- beta .* u[:,3]; 
end
prob = ODEProblem(lorenz!, u0, [0 1000], 1.0);


println("Ensemble")
# Remake simulatons varying a parameter (lambda), in parallel
function prob_func(prob, i, repeat);
    remake(prob, p = [lambdas[i]]);
end


# Solve the ODE in parallel calling the "remake" above
ensemble_prob = EnsembleProblem(prob, prob_func=prob_func);
sim = solve(ensemble_prob, Tsit5(), EnsembleThreads(), saveat = 900:0.25:1000, trajectories=len, maxiters = 1e5);


# Keep only the y coordinate
sim = copy(sim[:,1,:,:]);



######################################
########## COMPUTING ERRORS ##########
######################################

# Dictionary of cluster indices
cluster_indices = Dict("C1"=>[1,2,3], "C2"=>[4,5,6], "C3"=>[7,8,9,10])

# Call the functions from cluster_error_functions.jl
Err_matrix = compute_errors(cluster_indices, sim)



########################################
########## SAVING THE RESULTS ##########
########################################

num = rand(Int)
npzwrite("Output/simulation-$num.npy", Err_matrix)
