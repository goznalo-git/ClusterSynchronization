using DifferentialEquations
using MAT
using SparseArrays
using LinearAlgebra
using Statistics
using DelimitedFiles

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

length = 400 # Number of couplings (lambda) to calculate over 
lambdas = range(0, 8, length = length)


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
sim = solve(ensemble_prob, Tsit5(), EnsembleThreads(), saveat = 900:0.25:1000, trajectories=length, maxiters = 1e5);


# Keep only the y coordinate
sim = copy(sim[:,1,:,:]);



######################################
########## COMPUTING ERRORS ##########
######################################

# List of cluster indices
ind1 = [1,2,3];
ind2 = [4,5,6];
ind3 = [7,8,9,10];


# Compute the error corresponding to a single cluster
function ClusterError(soly);

    meany  = mean(soly, dims=1);
    diff = (soly .- meany).^2;
    sqroot = sqrt.(mean(diff, dims=1));
    Err  = mean(sqroot, dims=2)[:,1];
    
    return Err;
end


# Function computing the cluster error and total error of a single simulation
function clusynch(lam);

    # Error Synch Cluster 1
    listsol1 = [sim[x,:,lam] for x in ind1];
    soly1 = mapreduce(permutedims, vcat, listsol1[:][:]); #vector to matrix
    Err1 = ClusterError(soly1);

    # Error Synch Cluster 2
    listsol2 = [sim[x,:,lam] for x in ind2];
    soly2 = mapreduce(permutedims, vcat, listsol2[:][:]);
    Err2 = ClusterError(soly2);

    # Error Synch Cluster 3
    listsol3 = [sim[x,:,lam] for x in ind3];
    soly3 = mapreduce(permutedims, vcat, listsol3[:][:]);
    Err3 = ClusterError(soly3);
    
    # Error Synch Global
    solyT  = sim[:,:,lam];
    ErrT = ClusterError(solyT);
    
    return Err1, Err2, Err3, ErrT;
end


# Initialize cluster error vectors
x1 = Vector{Float64}()
x2 = Vector{Float64}()
x3 = Vector{Float64}()
xT = Vector{Float64}()

# For each value of lambda, append the cluster errors
for lam in 1:length;

    Err1, Err2, Err3, ErrT = clusynch(lam)
    append!(x1,Err1)
    append!(x2,Err2)
    append!(x3,Err3)
    append!(xT,ErrT)
end



########################################
########## SAVING THE RESULTS ##########
########################################

Synch = [x1; x2; x3; xT];

# Each different simulation is indexed by a random number
num = rand(Int)
open("Output/synch_error-$num.txt", "w") do file
    writedlm(file, Synch)
end
