"""
`load_model`

Load a model from the library.

**Arguments**
- `model_name`: Name of a model. Currently implemented models are:
    - `amplitude-phase`
    - `goodwin`
    - `goodwin-general` (kwargs: `n_equations`)
    - `van-der-pol`

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
        solver_parameters = (saveat=0.01, reltol=1e-9, abstol=1e-9,)

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
        
        # Parameters
        K, n, a1, a2, a3, d1, d2, d3 = p
        
        # Equations
        du[1] = a1 * K^n / (K^n + z^n) - d1 * x
        du[2] = a2 * x - d2 * y
        du[3] = a3 * y - d3 * z

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
        solver_parameters = (saveat=0.01, reltol=1e-9, abstol=1e-9,)

    else

        msg = "Problem type `$problem_type` not implemented! Use `ode`."
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
        solver_parameters = (saveat=0.01, reltol=1e-9, abstol=1e-9,)

    else

        msg = "Problem type `$problem_type` not implemented! Use `ode`."
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
        solver_parameters = (saveat=0.01, reltol=1e-9, abstol=1e-9,)

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
