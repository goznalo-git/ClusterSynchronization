
## This file contains the ODEproblem functions corresponding to Class III dynamical systems to be loaded in the ClusterSynchErrors scripts.


# Rossler ("official")
function rossler!(du, u, lambda, t, alpha=0.2, beta=0.2, gamma=9)
    coupling = similar(x0);
    mul!(coupling, sparseL, u[:,1]);
    du[:,1] = - u[:,2] .- u[:,3] .- lambda.*coupling;
    du[:,2] = u[:,1] .+ alpha * u[:,2];
    du[:,3] = beta .+ (u[:,1] .- gamma) .* u[:,3]; 
end
prob = ODEProblem(rossler!, u0, [0 10000], lam);


# Lorenz (beta=2, 3->3)
function lorenz!(du, u, lambda, t, sigma=10, rho=28, beta=2)
    coupling = similar(z0);
    mul!(coupling, sparseL, u[:,3]);
    du[:,1] = sigma * (u[:,2] .- u[:,1]);
    du[:,2] = u[:,1] .* (rho .- u[:,3]) .- u[:,2]; 
    du[:,3] = u[:,1] .* u[:,2] .- beta .* u[:,3] .- lambda .* coupling;
end
prob = ODEProblem(lorenz!, u0, [0 10000], lam);


# Lorenz (beta=2, 2->1)
function lorenz!(du, u, lambda, t, sigma=10, rho=28, beta=2)
    coupling = similar(y0);
    mul!(coupling, sparseL, u[:,2]);
    du[:,1] = sigma * (u[:,2] .- u[:,1]) .- lambda .* coupling;
    du[:,2] = u[:,1] .* (rho .- u[:,3]) .- u[:,2]; 
    du[:,3] = u[:,1] .* u[:,2] .- beta .* u[:,3];
end
prob = ODEProblem(lorenz!, u0, [0 10000], lam);


# Lorenz (beta=8/3, 2->1)
function lorenz!(du, u, lambda, t, sigma=10, rho=28, beta=8/3)
    coupling = similar(y0);
    mul!(coupling, sparseL, u[:,2]);
    du[:,1] = sigma * (u[:,2] .- u[:,1]) .- lambda .* coupling;
    du[:,2] = u[:,1] .* (rho .- u[:,3]) .- u[:,2]; 
    du[:,3] = u[:,1] .* u[:,2] .- beta .* u[:,3];
end
prob = ODEProblem(lorenz!, u0, [0 10000], lam);


# Chen
function chen!(du, u, lambda, t, a=35, c=28, beta=8/3)
    coupling = similar(z0);
    mul!(coupling, sparseL, u[:,3]);
    du[:,1] = a .* (u[:,2] .- u[:,1]);
    du[:,2] = (c .- a .- u[:,3]) .* u[:,1] .+ c .* u[:,2]; 
    du[:,3] = u[:,1] .* u[:,2] .- beta .* u[:,3] .- lambda .* coupling;
end
prob = ODEProblem(chen!, u0, [0 10000], lam);


# Duffing 
function duffing!(du, u, lambda, t, h=0.1, q=5.6, eta=1)
    coupling = similar(x0);
    mul!(coupling, sparseL, u[:,1]);
    du[:,1] = u[:,2];
    du[:,2] = - h .* u[:,2] .- u[:,1].^3 .+ q * sin(eta * t) .- lambda .* coupling;
end
prob = ODEProblem(duffing!, u0, [0 10000], lam);
