using DifferentialEquations
using MAT
using SparseArrays
using LinearAlgebra
using Statistics
using DelimitedFiles


# Load the appropriate adjacency matrix and convert it to Laplacian
#file1 = matread("Input/A_PowerGrid.mat");
#Aij = file1["A"];
#degree_vector = sum(Aij, dims=2);
#Lij = diagm(vec(degree_vector)) - Aij;

# or load directly the Laplacian matrix
file1 = matread("Input/L_SmallWeighted.mat");
Lij = file1["L"];

#Sparsify the Laplacian matrix
sparseL = sparse(Lij);


# Load the Rossler/Lorenz initial condition (within the attractor)
file2 = matread("Input/initc_Rossler.mat");
initc   = file2["initc"];

# In case of the 10% heterogeneity in variables a, b, c; 
# the values per node were already generated and saved in the corresponding files.
#het_b = vec(readdlm("Input/het_b.txt"));  
# If this is used, then the rossler!() function must be adapted accordingly.

nth = Threads.nthreads()
println(nth);


# Parameters and initial conditions
N = size(sparseL, 1);
x0 = initc[1] .+ 0.01*rand(N);
y0 = initc[2] .+ 0.01*rand(N);
z0 = initc[3] .+ 0.01*rand(N);
u0 = [x0 y0 z0];

length=200
lambdas = range(0, 0.2, length = length)


# ODE 
function rossler!(du, u, lambda, t, a=0.1, b=0.1, c=18)
    coupling = similar(y0);
    mul!(coupling, sparseL, u[:,2]);
    du[:,1] = - u[:,2] .- u[:,3];
    du[:,2] = u[:,1] .+ a * u[:,2] .- lambda.*coupling;
    du[:,3] = b .- c * u[:,3] .+ u[:,3] .* u[:,1]; 
end
prob = ODEProblem(rossler!,u0,[0 1000],1.0);


println("Ensemble")
# Remake simulatons varying a parameter: lambda
function prob_func(prob, i, repeat);
    remake(prob, p = [lambdas[i]]);
end

# Solve the ODE calling the "remake"
ensemble_prob = EnsembleProblem(prob, prob_func=prob_func);
sim = solve(ensemble_prob, Tsit5(), EnsembleThreads(), saveat = 900:0.25:1000, trajectories=length, maxiters = 1e5);

# Keep the y coordinate
sim = copy(sim[:,2,:,:]);


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


# For each value of lambda, compute the cluster errors
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


for lam in 1:length;

    Err1, Err2, Err3, ErrT = clusynch(lam)
    append!(x1,Err1)
    append!(x2,Err2)
    append!(x3,Err3)
    append!(xT,ErrT)
end

Synch = [x1; x2; x3; xT];

num = rand(Int)
open("Output-SmallWeighted/synch_error-$num.txt", "w") do file
    writedlm(file, Synch)
end
