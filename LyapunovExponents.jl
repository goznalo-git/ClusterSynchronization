#! /opt/julia/bin/julia 

using MAT
using LinearAlgebra
using DelimitedFiles

include("Lyapunov_functions.jl")

# Load the Lorenz initial condition (within the attractor)
file2 = matread("Input/initc_Rossler.mat");
initc = file2["initc"];
# Parameters
a   = 0.1;
b   = 0.1;
c   = 18;

# Range of coupling
nu  = 0:0.05:12;

# Time
t0      = 0;           # initial t
T       = 2000;        # Tmax
N       = 2000000;     # Number of iterations

t       = LinRange(t0, T, N);         # vector of time
h       = (T - t0)/N;  # delta t (spacing)

frec    = 50;    # frequency at which we measure/calculate the Lyapunov exponent CAREFUL WITH THIS!

Tsample = 15000;       # Late times used for the average


## Run the program
explya = LyapExp(nu, initc, t, h, N, frec, Tsample)


open("Rosslertest.txt", "w") do file
    writedlm(file, explya)
end
