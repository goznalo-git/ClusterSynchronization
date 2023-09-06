#! /opt/bin/julia 

using MAT
using LinearAlgebra


"""
    RK_step(f, params, ti, p, h)

Compute a single Runge-Kutta step with the function f, parameters params, 
at time ti with values given by the vector p and a spacing h.
"""
function RK_step(f, params, ti, p, h)
   
    k1 = f(ti, p, params);
    k2 = f(ti + 0.5 * h, p + 0.5 * h * k1, params);
    k3 = f(ti + 0.5 * h, p + 0.5 * h * k2, params);
    k4 = f(ti + h      , p + h * k3, params);
    
    return p + (h/6) * (k1 + 2 * k2 + 2 * k3 + k4);
end



"""
    compute_RK4(f, params, p0, t, h, N, frec)

Compute a Runge-Kutta step with the function f, with parameters params, during time t, with timestep h
starting from vector p0, over N iterations, and reset after frec iterations.
"""
function compute_RK4(f, params, p0, t, h, N, frec)
    P       = zeros(N,6);                 # Initialize the eqs (6, deviation and standard) vector P 
    P[1,:]  = p0;                         # Initial condition
    tl      = zeros(Int(N/frec - 1));  # vector of time at which the measures were taken.
    enorm   = zeros(Int(N/frec - 1));                # initialization of the error norm's vector
    
    # Runge-Kutta
    index = 1;
    for it = 1:N-1
        
        # Take the new p, t values and perform a single step
        p  = P[it,:];
        ti = t[it];
        P[it + 1, :] = RK_step(f, params, ti, p, h)

        # Measure of the error vector norm every frec iterations
        if mod(it,frec) == 0
            enorm[index] = log(norm(P[it, 1:3]))/(frec * h); # log of error vector norm divided by elapsed time 
            P[it + 1, 1:3] = normalize(P[it + 1, 1:3]);      # Normalize the error vector for the next iteration
            
            tl[index] = it * h;                              
            
            index += 1;
        end
    end
    
    return P, enorm, tl
end


"""
    LyapExp(nulist, f, parameters, initc, t, h, N, frec, Tsample)

Compute the Lyapunov exponent of the dynamical system f(x) with parameters,
couplings nulist, from an initial condition initc, during N iterations with h timestep, 
calculated at a frequency frec, and averaged over the last Tsample values.
"""
function LyapExp(nulist, f, parameters, initc, t, h, N, frec, Tsample)
    
    println(parameters[1])

    explya = zeros(length(nulist))
    for (index, nu) in enumerate(nulist)
        
        # initialization  p0(1:3,1) random. p0(4:6,1) initc
        p0          = zeros(6);    
        p0[4:6, 1]  = ones(3) .* initc';  # Decoupled system
        p0[1:3, 1]  = normalize(rand(3)); # Error vector
        
        # E.g. the a,b,c of Rossler or the sigma, rho, beta of Lorenz
        params = (parameters[1], parameters[2], parameters[3], nu)
        
        P, enorm, tl = compute_RK4(f, params, p0, t, h, N, frec);
        
        explya[index]   = sum(enorm[end - Tsample:end]) / tl[end];

    end
    
    return explya
end