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

# Load the appropriate adjacency matrix and convert it to Laplacian
file1 = matread("Input/A_PowerGrid.mat");
Aij = file1["A"];
degree_vector = sum(Aij, dims=2);
Lij = diagm(vec(degree_vector)) - Aij;

#Sparsify the Laplacian matrix
sparseL = sparse(Lij);

# Load the Rossler initial condition (within the attractor)
file2 = matread("Input/initc_Rossler.mat");
initc   = file2["initc"];


##########################################
########## ESTABLISH PARAMETERS ##########
##########################################


# Parameters and initial conditions
N = size(sparseL, 1);
x0 = initc[1] .+ 0.01*rand(N);
y0 = initc[2] .+ 0.01*rand(N);
z0 = initc[3] .+ 0.01*rand(N);
u0 = [x0 y0 z0];

length = 200 # Number of couplings (lambda) to calculate over 
lambdas = range(0, 0.2, length = length)


# Check the number of threads available for the parallelization
nth = Threads.nthreads()
println(nth);


##########################################
########## EQUATIONS AND SOLVER ##########
##########################################

# Establish ODE 
function rossler!(du, u, lambda, t, a=0.1, b=0.1, c=18)
    coupling = similar(y0);
    mul!(coupling, sparseL, u[:,2]);
    du[:,1] = - u[:,2] .- u[:,3];
    du[:,2] = u[:,1] .+ a * u[:,2] .- lambda.*coupling;
    du[:,3] = b .- c * u[:,3] .+ u[:,3] .* u[:,1]; 
end
prob = ODEProblem(rossler!,u0,[0 1000],1.0);


println("Ensemble")
# Remake simulatons varying a parameter (lambda), in parallel
function prob_func(prob, i, repeat);
    remake(prob, p = [lambdas[i]]);
end


# Solve the ODE in parallel calling the "remake" above
ensemble_prob = EnsembleProblem(prob, prob_func=prob_func);
sim = solve(ensemble_prob, Tsit5(), EnsembleThreads(), saveat = 900:0.25:1000, trajectories=length, maxiters = 1e5);


# Keep only the y coordinate
sim = copy(sim[:,2,:,:]);



######################################
########## COMPUTING ERRORS ##########
######################################

# List of cluster indices
cluster_indices = Dict("C1"=>[1,2,3], "C2"=>[4,5,6], "C3"=>[7,8,9,10])

# Call the functions from cluster_error_functions.jl
Err_matrix = compute_errors(cluster_indices, sim)



########################################
########## SAVING THE RESULTS ##########
########################################

num = rand(Int)
npzwrite("Output/error_test.npy", Err_matrix)
