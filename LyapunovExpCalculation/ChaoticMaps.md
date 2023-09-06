# Different Chaotic Maps for the Lyapunov exponents calculation 

These must be inserted in the `LyapunovExponents.jl` script, along with the specified parameters.

Only a partial list (the ones used in some of our works) are available here.

1. Rossler

    - Parameters sets
    ```julia
    a = 0.2;
    b = 0.2;
    c = 9;
    
    parameters = (a, b, c)
    ```
    
    - Coupling 1 -> 1
    ```julia
    function f(t, p, params)

        dx, dy, dz, x, y, z = p
        a, b, c, nu = params

        return [- dy - dz - nu * dx;  dx + a * dy;  z * dx + (x - c) * dz;  # error
                            - y - z;    x + a * y;       b + (x - c) * z];  # decoupled
    end
    ```


2. Lorenz

    - Parameters sets
    ```julia
    sigma = 10;
    rho   = 28;
    beta  = 2;
    
    sigma = 10;
    rho   = 28;
    beta  = 8/3;
    
    parameters = (sigma, rho, beta)
    ```

    - Coupling 2 -> 1
    ```julia
    function f(t, p, params)

        dx, dy, dz, x, y, z = p
        sigma, rho, beta, nu = params

        return [sigma * (dy - dx) - nu * dy; (rho - z) * dx - dy - x * dz; y * dx + x * dy - beta * dz;  # error
                            sigma * (y - x);            x * (rho - z) - y;           x * y - beta * z];  # decoupled 
    end
    ```

3. Chen

    - Parameters sets
    ```julia
    a = 35;
    c   = 28;
    beta  = 8/3;
    
    parameters = (a, c, beta)
    ```

    - Coupling 3 -> 3
    ```julia
    function f(t, p, params)

        dx, dy, dz, x, y, z = p
        sigma, rho, beta, nu = params

        return [a * (dy - dx);  (c - a - z) * dx - dz * x - c * dy; y * dx + x * dy - beta * dz - nu * dz;  # error
                  a * (y - x);             (c - a - z) * x - c * y;                     x * y - beta * z];  # decoupled 
    end
    ```

4. Chua

    - Parameter sets
    ```julia
    alpha = 10
    beta = 14.87
    gamma = 0
    a = -1.27
    b = -0.68
    
    parameters = (alpha, beta, gamma, a, b)
    ```

    - Coupling 3 -> 3

    ```
    function f(t, p, params)
        
        dx, dy, dz, x, y, z = p
        alpha, beta, gamma, a, b = params
        
        if x > 1
            
            fx = - b * x - a + b
            dfx = - b

        elif x < -1
            
            fx = - b * x + a - b
            dfx = - b

        else
            
            fx = - a * x
            dfx = - a

        end

        return [alpha * (dy - dx + dfx);  dx - dy + dz; - beta * dy - gamma * dz - beta * dz - nu * dz;  # error
                  alpha * (y - x);             (c - a - z) * x - c * y;                     x * y - beta * z];  # decoupled 
    end
    ```
