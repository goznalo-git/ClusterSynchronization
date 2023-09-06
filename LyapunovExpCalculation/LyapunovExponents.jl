#! /opt/julia/bin/julia 

using MAT
using LinearAlgebra
using DelimitedFiles

include("Lyapunov_functions.jl")


println("Starting the calculation")

# Load the Lorenz initial condition (within the attractor)
file2 = matread("Input/initc_Rossler.mat");
initc = file2["initc"];

# Parameters
sigma   = 10;
rho     = 28;
beta    = 2;

parameters = (sigma, rho, beta)

# Range of coupling
nulist  = 0:0.5:25;

# Time
t0      = 0;           # initial t
T       = 2000;        # Tmax
N       = 2000000;     # Number of iterations

t       = LinRange(t0, T, N);         # vector of time
h       = (T - t0)/N;  # delta t (spacing)

frec    = 50;    # frequency at which we measure/calculate the Lyapunov exponent CAREFUL WITH THIS!

Tsample = 15000;       # Late times used for the average



###########################################
######## Taken from ChaoticMaps.md ########
###########################################

function f(t, p, params)

    dx, dy, dz, x, y, z = p
    sigma, rho, beta, nu = params

    return [sigma * (dy - dx) - nu * dy; (rho - z) * dx - dy - x * dz; y * dx + x * dy - beta * dz;  # error
                        sigma * (y - x);            x * (rho - z) - y;           x * y - beta * z];  # decoupled 
end

###########################################
###########################################

## Run the program
explya = LyapExp(nulist, f, parameters, initc, t, h, N, frec, Tsample)


open("Lorenz2-1.txt", "w") do file
    writedlm(file, explya)
end
