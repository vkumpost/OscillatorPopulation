"""
`load_model`

Load a model from the library.

**Arguments**
- `model_name`: Name of a model. Currently implemented models are:
    - `amplitude-phase`
    - `goodwin`
    - `goodwin-general` (kwargs: `n_equations`; default value is 3)
    - `van-der-pol`
    - `kim-forger`
    - `kim-forger-full`
    - `hasty`

**Optional Arguments**
- `problem_type`: "ode", "sde", or "jump".
"""
function load_model(model_name::String, problem_type::String="ode"; kwargs...)

    if model_name == "amplitude-phase"
        model = _load_amplitude_phase(problem_type)
    elseif model_name == "goodwin"
        model = _load_goodwin(problem_type)
    elseif model_name == "goodwin-general"
        model = _load_goodwin_general(problem_type; kwargs...)
    elseif model_name == "van-der-pol"
        model = _load_van_der_pol(problem_type)
    elseif model_name == "kim-forger"
        model = _load_kim_forger(problem_type)
    elseif model_name == "kim-forger-full"
        model = _load_kim_forger_full(problem_type)
    elseif model_name == "hasty"
        model = _load_hasty(problem_type)
    else
        err = OscillatorPopulationError("`$model_name` is not in the model library!")
        throw(err)
    end

    return model

end


"""
`_load_amplitude_phase`

Load the amplitude-phase model.

**References**
- Isao T. Tokuda, et al. "Conceptual models of entrainment, jet lag, and
    seasonality." Frontiers in Physiology 11 (2020): 334.
"""
function _load_amplitude_phase(problem_type)

    function model!(du, u, p, t)

        # Variables
        x, y = u
        
        # Parameters
        λ, A, T, I = p
        
        # Equations
        r = sqrt(x^2 + y^2)
        ω = 2π/T
        du[1] = -λ * x * (r - A) - ω * y + I
        du[2] = -λ * y * (r - A) + ω * x

    end

    function noise!(du, u, p, t)
    
        # Parameters
        σ = p[end]
    
        # Equations
        du[1] = σ
        du[2] = σ

    end

    # Common variables for all model types
    tspan = (0.0, 10.0)
    variable_names = ["x", "y"]
    parameter_names = ["λ", "A", "T", "I"]
    p = [0.5, 2.0, 1.0, 0.0]
    u0 = [0.1, 0.1]
    input = (Matrix{Float64}(undef, 0, 0), "")
    output = sol -> Matrix(sol[:, :]')
        
    if problem_type == "ode"
        
        # Create an ODE problem
        problem = ODEProblem(model!, u0, tspan, p)
    
        # Set solver parameters
        solver_algorithm = DP5()
        solver_parameters = (saveat=0.01,)

    elseif problem_type == "sde"

        # Extend parameters by the noise intensity
        push!(parameter_names, "σ")
        push!(p, 0.5)

        # Create an SDE problem
        problem = SDEProblem(model!, noise!, u0, tspan, p)
                
        # Set solver parameters
        solver_algorithm = SOSRI()
        solver_parameters = (saveat=0.01,)

    else

        msg = "Problem type `$problem_type` not implemented! Use `ode` or `sde`."
        err = OscillatorPopulationError(msg)
        throw(err)

    end

    model = Model(variable_names, parameter_names, problem, solver_algorithm,
        solver_parameters, input, output)
        
    return model

end


"""
`_load_goodwin`

Load the Goodwin model.

**References**
- Didier Gonze and Peter Ruoff. "The Goodwin oscillator and its legacy." Acta
    Biotheoretica 69.4 (2021): 857-874.
"""
function _load_goodwin(problem_type)

    function model!(du, u, p, t)

        # Variables
        x, y, z = u
        x = max(0, x)
        y = max(0, y)
        z = max(0, z)
        
        # Parameters
        K, n, a1, a2, a3, d1, d2, d3 = p
        
        # Equations
        du[1] = a1 * K^n / (K^n + z^n) - d1 * x
        du[2] = a2 * x - d2 * y
        du[3] = a3 * y - d3 * z

    end

    function noise!(du, u, p, t)

        # Variables
        x, y, z = u
        x = max(0, x)
        y = max(0, y)
        z = max(0, z)

        # Parameters
        K, n, a1, a2, a3, d1, d2, d3, σ = p

        # Equations
        du[1, 1] = σ * sqrt(a1 * K^n / (K^n + z^n))
        du[1, 2] = σ * -sqrt(d1 * x)
        du[2, 3] = σ * sqrt(a2 * x)
        du[2, 4] = σ * -sqrt(d2 * y)
        du[3, 5] = σ * sqrt(a3 * y)
        du[3, 6] = σ * -sqrt(d3 * z)

    end

    # Common variables for all model types
    tspan = (0.0, 100.0)
    variable_names = ["x", "y", "z"]
    parameter_names = ["K", "n", "a1", "a2", "a3", "d1", "d2", "d3"]
    p = [1.0, 10.0, 5.0, 5.0, 5.0, 0.5, 0.5, 0.5]
    u0 = [0.1, 0.1, 0.1]
    input = (Matrix{Float64}(undef, 0, 0), "")
    output = sol -> Matrix(sol[:, :]')
        
    if problem_type == "ode"
        
        # Create an ODE problem
        problem = ODEProblem(model!, u0, tspan, p)
    
        # Set solver parameters
        solver_algorithm = DP5()
        solver_parameters = (saveat=0.01,)

    elseif problem_type == "sde"

        # Create an SDE model
        push!(parameter_names, "σ")
        push!(p, 0.1)
        noise_rate_prototype = zeros(3, 6)
        problem = SDEProblem(model!, noise!, u0, tspan, p;
            noise_rate_prototype=noise_rate_prototype)

        # Set solver parameters
        solver_algorithm = EM()
        solver_parameters = (dt=0.0001, saveat=0.01,)

    else

        msg = "Problem type `$problem_type` not implemented! Use `ode` or `sde`."
        err = OscillatorPopulationError(msg)
        throw(err)

    end

    model = Model(variable_names, parameter_names, problem, solver_algorithm,
        solver_parameters, input, output)
        
    return model

end


"""
`_load_goodwin_general`

Load the generalized Goodwin model.

**References**
- Sanchez, Luis A. "Global asymptotic stability of the Goodwin system with
    repression." Nonlinear analysis: real world applications 10.4 (2009):
    2151-2156.
"""
function _load_goodwin_general(problem_type; n_equations=3)

    function model!(du, u, p, t)
        
        # Parameters
        K, n = p
        
        # Equations
        du[1] = K^n / (K^n + u[end]^n) - u[1]
        for i = 2:n_equations
            du[i] = u[i-1] - u[i]
        end

    end

    function noise!(du, u, p, t)

        # Parameters
        K, n, σ = p

        # Equations
        for i = 1:n_equations
            du[i] = σ
        end

    end

    # Common variables for all model types
    tspan = (0.0, 10.0)
    variable_names = ["x$i" for i in 1:n_equations]
    parameter_names = ["K", "n"]
    p = [0.1, 12.0]
    u0 = fill(0.1, n_equations)
    input = (Matrix{Float64}(undef, 0, 0), "")
    output = sol -> Matrix(sol[:, :]')
        
    if problem_type == "ode"
        
        # Create an ODE problem
        problem = ODEProblem(model!, u0, tspan, p)
    
        # Set solver parameters
        solver_algorithm = DP5()
        solver_parameters = (saveat=0.01,)

    elseif problem_type == "sde"

        # Create an SDE model
        push!(parameter_names, "σ")
        push!(p, 0.01)
        problem = SDEProblem(model!, noise!, u0, tspan, p)

        # Set solver parameters
        solver_algorithm = SOSRI()
        solver_parameters = (saveat=0.01,)

    else

        msg = "Problem type `$problem_type` not implemented! Use `ode` or `sde`."
        err = OscillatorPopulationError(msg)
        throw(err)

    end

    model = Model(variable_names, parameter_names, problem, solver_algorithm,
        solver_parameters, input, output)
        
    return model

end


"""
`_load_van_der_pol`

Load the Van der Pol model.

**References**
- Namiko Mitarai, Uri Alon, and Mogens H. Jensen. "Entrainment of noise-induced
    and limit cycle oscillators under weak noise." Chaos: An Interdisciplinary
    Journal of Nonlinear Science 23.2 (2013): 023125.
"""
function _load_van_der_pol(problem_type)

    function model!(du, u, p, t)

        # Variables
        x, y = u
        
        # Parameters
        B, d, I = p
        
        # Equations
        du[1] = y
        du[2] = -(B * x^2 - d) * y - x + I

    end

    function noise!(du, u, p, t)
    
        # Parameters
        σ = p[end]
    
        # Equations
        du[1] = σ
        du[2] = σ

    end

    # Common variables for all model types
    tspan = (0.0, 100.0)
    variable_names = ["x", "y"]
    parameter_names = ["B", "d", "I"]
    p = [10.0, 2.0, 0.0]
    u0 = [0.1, 0.1]
    input = (Matrix{Float64}(undef, 0, 0), "")
    output = sol -> Matrix(sol[:, :]')
        
    if problem_type == "ode"
        
        # Create an ODE problem
        problem = ODEProblem(model!, u0, tspan, p)
    
        # Set solver parameters
        solver_algorithm = DP5()
        solver_parameters = (saveat=0.01,)

    elseif problem_type == "sde"

        # Extend parameters by the noise intensity
        push!(parameter_names, "σ")
        push!(p, 0.1)

        # Create an SDE problem
        problem = SDEProblem(model!, noise!, u0, tspan, p)
                
        # Set solver parameters
        solver_algorithm = SOSRI()
        solver_parameters = (saveat=0.01,)

    else

        msg = "Problem type `$problem_type` not implemented! Use `ode` or `sde`."
        err = OscillatorPopulationError(msg)
        throw(err)

    end

    model = Model(variable_names, parameter_names, problem, solver_algorithm,
        solver_parameters, input, output)
        
    return model

end


"""
`_load_kim_forger`

Load the Kim-Forger model with simplified parameter set.

**References**
- Jae Kyoung Kim, Daniel B. Forger. "A mechanism for robust circadian
    timekeeping via stoichiometric balance." Molecular Systems Biology 8 (2012).
    https://doi.org/10.1038/msb.2012.62
"""
function _load_kim_forger(problem_type)

    function model!(du, u, p, t)

        # Variables
        x, y, z = u
        
        # Parameters
        A, I, τ = p
        
        # Equations
        du[1] = τ * max(0, 1 - z/A) - τ * x + I
        du[2] = τ * x - τ * y
        du[3] = τ * y - τ * z

    end

    function noise!(du, u, p, t)

        # Variables
        x, y, z = u
        x = max(0, x)
        y = max(0, y)
        z = max(0, z)

        # Parameters
        A, I, τ, σ = p

        # Equations
        du[1, 1] = σ * sqrt(τ * max(0, 1 - z/A))
        du[1, 2] = σ * sqrt(τ * x)
        du[1, 3] = σ * sqrt(I)
        du[2, 4] = σ * sqrt(τ * x)
        du[2, 5] = σ * sqrt(τ * y)
        du[3, 6] = σ * sqrt(τ * y)
        du[3, 7] = σ * sqrt(τ * z)

    end

    function jumps()

        ix, iy, iz = 1:3
        iA, iI, iτ, iΩ = 1:4

        rate1(u, p, t) = p[iτ]*p[iΩ]*max(0, 1 - u[iz] / (p[iΩ] * p[iA]))
        affect1!(integrator) = integrator.u[ix] += 1
        jump1 = ConstantRateJump(rate1, affect1!)
    
        rate2(u, p, t) = p[iτ]*u[ix]
        affect2!(integrator) = integrator.u[ix] -= 1
        jump2 = ConstantRateJump(rate2, affect2!)
    
        rate3(u, p, t) = p[iτ]*u[ix]
        affect3!(integrator) = integrator.u[iy] += 1
        jump3 = ConstantRateJump(rate3, affect3!)
    
        rate4(u, p, t) = p[iτ]*u[iy]
        affect4!(integrator) = integrator.u[iy] -= 1
        jump4 = ConstantRateJump(rate4, affect4!)
    
        rate5(u, p, t) = p[iτ]*u[iy]
        affect5!(integrator) = integrator.u[iz] += 1
        jump5 = ConstantRateJump(rate5, affect5!)
    
        rate6(u, p, t) = p[iτ]*u[iz]
        affect6!(integrator) = integrator.u[iz] -= 1
        jump6 = ConstantRateJump(rate6, affect6!)
    
        rate7(u, p, t) = p[iΩ]*p[iI]
        affect7!(integrator) = (integrator.u[iX] += 1)
        jump7 = ConstantRateJump(rate7, affect7!)
    
        return jump1, jump2, jump3, jump4, jump5, jump6, jump7

    end

    # Common variables for all model types
    tspan = (0.0, 100.0)
    variable_names = ["x", "y", "z"]
    parameter_names = ["A", "I", "τ"]
    p = [0.1, 0.0, 3.66]
    u0 = [0.1, 0.1, 0.1]
    input = (Matrix{Float64}(undef, 0, 0), "")
    output = sol -> Matrix(sol[:, :]')
        
    if problem_type == "ode"
        
        # Create an ODE problem
        problem = ODEProblem(model!, u0, tspan, p)
    
        # Set solver parameters
        solver_algorithm = DP5()
        solver_parameters = (saveat=0.01,)

    elseif problem_type == "sde"

        # Create an SDE model
        push!(parameter_names, "σ")
        push!(p, 0.1)
            
        noise_rate_prototype = zeros(3, 7)
        problem = SDEProblem(model!, noise!, u0, tspan, p;
            noise_rate_prototype=noise_rate_prototype
        )

        # Set solver parameters
        solver_algorithm = EM()
        solver_parameters = (dt=0.001, saveat=0.01,)

    elseif problem_type == "jump"

        # Create a Jump model
        push!(parameter_names, "Ω")
        push!(p, 1000.0)
        u0 .*= p[end]
        discrete_problem = DiscreteProblem(u0, tspan, p)
        problem = JumpProblem(discrete_problem, Direct(), jumps()...;
            save_positions=(false, false))
            
        # Set solver parameters
        solver_algorithm = SSAStepper()
        solver_parameters = (saveat=0.01,)

    else

        msg = "Problem type `$problem_type` not implemented! Use `ode`, `sde`, or `jump`."
        err = OscillatorPopulationError(msg)
        throw(err)

    end

    model = Model(variable_names, parameter_names, problem, solver_algorithm,
        solver_parameters, input, output)
        
    return model


end


"""
`_load_kim_forger`

Load the Kim-Forger model.

**References**
- Jae Kyoung Kim, Daniel B. Forger. "A mechanism for robust circadian
    timekeeping via stoichiometric balance." Molecular Systems Biology 8 (2012).
    https://doi.org/10.1038/msb.2012.62
"""
function _load_kim_forger_full(problem_type)

    kfr = (R, A, K) -> (A - R - K + sqrt((A - R - K)^2 + 4*A*K)) / (2*A)

    function model!(du, u, p, t)

        # Variables
        M, P, R = u

        # Parameters
        A, K, vM, vP, vR, dM, dP, dR, I = p

        # Equations
        du[1] = dM = vM*kfr(R, A, K) - dM*M + I
        du[2] = dP = vP*M - dP*P
        du[3] = dR = vR*P - dR*R

    end

    function noise!(du, u, p, t)

        # Variables
        M, P, R = u
        M = max(0, M)
        P = max(0, P)
        R = max(0, R)

        # Parameters
        A, K, vM, vP, vR, dM, dP, dR, I, σ = p

        # Equations
        du[1, 1] = σ * sqrt(vM*kfr(R, A, K))
        du[1, 2] = σ * sqrt(dM*M)
        du[1, 3] = σ * sqrt(I)
        du[2, 4] = σ * sqrt(vP*M)
        du[2, 5] = σ * sqrt(dP*P)
        du[3, 6] = σ * sqrt(vR*P)
        du[3, 7] = σ * sqrt(dR*R)

    end

    # Common variables for all model types
    tspan = (0.0, 100.0)
    variable_names = ["M", "P", "R"]
    parameter_names = ["A", "K", "vM", "vP", "vR", "dM", "dP", "dR", "I"]
    p = [9.0, 0.0, 1.0, 1.0, 1.0, 0.16, 0.16, 0.16, 0.05]
    u0 = [0.1, 0.1, 0.1]
    input = (Matrix{Float64}(undef, 0, 0), "")
    output = sol -> Matrix(sol[:, :]')
        
    if problem_type == "ode"
        
        # Create an ODE problem
        problem = ODEProblem(model!, u0, tspan, p)
    
        # Set solver parameters
        solver_algorithm = DP5()
        solver_parameters = (saveat=0.01,)

    elseif problem_type == "sde"

        # Create an SDE model
        push!(parameter_names, "σ")
        push!(p, 0.04)
            
        noise_rate_prototype = zeros(3, 7)
        problem = SDEProblem(model!, noise!, u0, tspan, p;
            noise_rate_prototype=noise_rate_prototype
        )

        # Set solver parameters
        solver_algorithm = EM()
        solver_parameters = (dt=0.001, saveat=0.01,)

    else

        msg = "Problem type `$problem_type` not implemented! Use `ode` or `sde`."
        err = OscillatorPopulationError(msg)
        throw(err)

    end

    model = Model(variable_names, parameter_names, problem, solver_algorithm,
        solver_parameters, input, output)
        
    return model


end


"""
`_load_hasty`

Load the Hasty model.

**References**
- Hasty J, Dolnik M, Rottschäfer V, Collins JJ. "Synthetic gene network for
    entraining and amplifying cellular oscillations." Physical Review Letters 88
    (2002). https://doi.org/10.1103/PhysRevLett.88.148101
"""
function _load_hasty(problem_type)

    function model!(du, u, p, t)

        # Variables
        x, y = u

        # Parameters
        α, β, ay, dx, dy, I = p

        # Equations
        F = (1 + x^2 + α*β*x^4) / ((1 + x^2 + β*x^4)*(1 + y^4))
        du[1] = dx = F -  dx*x + I
        du[2] = dy = ay * F - dy*y

    end

    function noise!(du, u, p, t)

        # Variables
        x, y = u
        
        # Parameters
        α, β, ay, dx, dy, I, σ = p

        # Equations
        du[1] = σ
        du[2] = σ

    end

    # Common variables for all model types
    tspan = (0.0, 300.0)
    variable_names = ["x", "y"]
    parameter_names = ["α", "β", "ay", "dx", "dy", "I"]
    p = [11.0, 2.0, 0.2, 0.2, 0.012, 0.0]
    u0 = [0.1, 0.1]
    input = (Matrix{Float64}(undef, 0, 0), "")
    output = sol -> Matrix(sol[:, :]')
        
    if problem_type == "ode"
        
        # Create an ODE problem
        problem = ODEProblem(model!, u0, tspan, p)
    
        # Set solver parameters
        solver_algorithm = DP5()
        solver_parameters = (saveat=0.01,)

    elseif problem_type == "sde"

        # Extend parameters by the noise intensity
        push!(parameter_names, "σ")
        push!(p, 0.01)

        # Create an SDE problem
        problem = SDEProblem(model!, noise!, u0, tspan, p)
                
        # Set solver parameters
        solver_algorithm = SOSRI()
        solver_parameters = (saveat=0.01,)

    else

        msg = "Problem type `$problem_type` not implemented! Use `ode` or `sde`."
        err = OscillatorPopulationError(msg)
        throw(err)

    end

    model = Model(variable_names, parameter_names, problem, solver_algorithm,
        solver_parameters, input, output)
        
    return model

end
