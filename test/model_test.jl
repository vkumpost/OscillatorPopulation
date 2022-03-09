"""
Create a dummy model for testing.
"""
function create_dummy(problem_type)

    function model!(du, u, p, t)
        du[1] = p[1]
        du[2] = p[2]
    end

    function noise!(du, u, p, t)
        du[1] = p[3]
        du[2] = p[3]
    end

    function jumps()

        rate1(u, p, t) = p[1]
        affect1!(integrator) = integrator.u[1] += 1
        jump1 = ConstantRateJump(rate1, affect1!)

        rate2(u, p, t) = p[2]
        affect2!(integrator) = integrator.u[2] += 1
        jump2 = ConstantRateJump(rate2, affect2!)

        return jump1, jump2

    end

    variable_names = ["x", "y"]
    parameter_names = ["a", "b", "c"]
    u0 = [1.0, 0.0]
    tspan = (0.0, 4.0)
    p = [1.0, 2.0, 0.0]

    if problem_type == "ode"
        problem = ODEProblem(model!, u0, tspan, p)
        solver_algorithm = Euler()
        solver_parameters = (dt=1.0,)

    elseif problem_type == "sde"
        problem = SDEProblem(model!, noise!, u0, tspan, p)
        solver_algorithm = EM()
        solver_parameters = (dt=1.0,)

    elseif problem_type == "jump"
        discrete_problem = DiscreteProblem(u0, tspan, p)
        problem = JumpProblem(discrete_problem, Direct(), jumps()...;
            save_positions=(false, false))
        solver_algorithm = SSAStepper()
        solver_parameters = (saveat=1.0,)

    end

    input = (Matrix{Float64}(undef, 0, 0), "")
    output = sol -> Matrix(sol[:, :]')
    
    model = OscillatorPopulation.Model(variable_names, parameter_names, problem,
        solver_algorithm, solver_parameters, input, output)

    return model

end


@testset "create_dummy" begin
    
    # ode model
    model = create_dummy("ode")
    @test model.variable_names == ["x", "y"]
    @test model.parameter_names == ["a", "b", "c"]
    @test model.problem isa ODEProblem
    @test model.problem.u0 == [1.0, 0.0]
    @test model.problem.tspan == (0.0, 4.0)
    @test model.problem.p == [1.0, 2.0, 0.0]
    @test model.solver_algorithm == Euler()
    @test model.solver_parameters == (dt=1.0,)
    @test isempty(model.input[1])
    @test model.input[1] isa Matrix
    @test isempty(model.input[2])
    @test model.input[2] isa String
    @test model.output isa Function

    # sde model
    model = create_dummy("sde")
    @test model.variable_names == ["x", "y"]
    @test model.parameter_names == ["a", "b", "c"]
    @test model.problem isa SDEProblem
    @test model.problem.u0 == [1.0, 0.0]
    @test model.problem.tspan == (0.0, 4.0)
    @test model.problem.p == [1.0, 2.0, 0.0]
    @test model.solver_algorithm == EM()
    @test model.solver_parameters == (dt=1.0,)
    @test isempty(model.input[1])
    @test model.input[1] isa Matrix
    @test isempty(model.input[2])
    @test model.input[2] isa String
    @test model.output isa Function

    # jump model
    model = create_dummy("jump")
    @test model.variable_names == ["x", "y"]
    @test model.parameter_names == ["a", "b", "c"]
    @test model.problem isa JumpProblem
    @test model.problem.prob.u0 == [1.0, 0.0]
    @test model.problem.prob.tspan == (0.0, 4.0)
    @test model.problem.prob.p == [1.0, 2.0, 0.0]
    @test model.solver_algorithm == SSAStepper()
    @test model.solver_parameters == (saveat=1.0,)
    @test isempty(model.input[1])
    @test model.input[1] isa Matrix
    @test isempty(model.input[2])
    @test model.input[2] isa String
    @test model.output isa Function

end


@testset "_get_problem_property" begin

    model = create_dummy("ode")
    u0 = OscillatorPopulation._get_problem_property(model, "u0")
    p = OscillatorPopulation._get_problem_property(model, "p")
    tspan = OscillatorPopulation._get_problem_property(model, "tspan")
    @test u0 == [1.0, 0.0]
    @test p == [1.0, 2.0, 0.0]
    @test tspan == (0.0, 4.0)

    model = create_dummy("jump")
    u0 = OscillatorPopulation._get_problem_property(model, "u0")
    p = OscillatorPopulation._get_problem_property(model, "p")
    tspan = OscillatorPopulation._get_problem_property(model, "tspan")
    @test u0 == [1.0, 0.0]
    @test p == [1.0, 2.0, 0.0]
    @test tspan == (0.0, 4.0)

end

@testset "_set_problem_property!" begin

    model = create_dummy("ode")
    OscillatorPopulation._set_problem_property!(model, u0=[0.0, 0.0], p=[0.0, 0.0, 2.0])
    @test model.problem.u0 == [0.0, 0.0]
    @test model.problem.p == [0.0, 0.0, 2.0]
    @test model.problem.tspan == (0.0, 4.0)

    model = create_dummy("sde")
    OscillatorPopulation._set_problem_property!(model, p=[0.0, 0.0, 2.0])
    @test model.problem.u0 == [1.0, 0.0]
    @test model.problem.p == [0.0, 0.0, 2.0]
    @test model.problem.tspan == (0.0, 4.0)

    model = create_dummy("jump")
    OscillatorPopulation._set_problem_property!(model, tspan=(0.0, 12.3), u0=[12.0, -5.0])
    @test model.problem.prob.u0 == [12.0, -5.0]
    @test model.problem.prob.p == [1.0, 2.0, 0.0]
    @test model.problem.prob.tspan == (0.0, 12.3)

end

@testset "set_initial_conditions!" begin
    
    model = create_dummy("ode")
    set_initial_conditions!(model, [3.1, -1.8])
    @test model.problem.u0 == [3.1, -1.8]

    model = create_dummy("sde")
    set_initial_conditions!(model, [3.1, -1.8])
    @test model.problem.u0 == [3.1, -1.8]

    model = create_dummy("jump")
    set_initial_conditions!(model, [3.1, -1.8])
    @test model.problem.prob.u0 == [3.1, -1.8]

end

@testset "set_timespan!" begin

    model = create_dummy("ode")
    set_timespan!(model, (-1.1, 49.2))
    @test model.problem.tspan == (-1.1, 49.2)

    model = create_dummy("sde")
    set_timespan!(model, (-1.1, 49.2))
    @test model.problem.tspan == (-1.1, 49.2)

    model = create_dummy("jump")
    set_timespan!(model, (-1.1, 49.2))
    @test model.problem.prob.tspan == (-1.1, 49.2)

end

@testset "set_output!" begin

    model = create_dummy("ode")
    fun = sol -> sol.t
    set_output!(model, fun)
    @test model.output == fun

end

@testset "set_solver!" begin
    
    # Change solver and add a solver parameter
    model = create_dummy("ode")
    set_solver!(model, DP5(), saveat=5.0)
    @test model.solver_algorithm == DP5()  # changed
    @test model.solver_parameters.dt == 1.0  # not changed
    @test model.solver_parameters.saveat == 5.0  # added

    # Keep solver and change an existing solver parameter
    model = create_dummy("ode")
    set_solver!(model, dt=0.01)
    @test model.solver_algorithm == Euler()  # not change
    @test model.solver_parameters.dt == 0.01  # changed

    # Replace solver parameters with new ones
    model = create_dummy("ode")
    set_solver!(model, saveat=0.1, merge_kwargs=false)
    @test model.solver_algorithm == Euler()  # not changed
    @test model.solver_parameters.saveat == 0.1  # added
    @test !(:dt in keys(model.solver_parameters))  # removed

    # Check that callback is not removed
    model = create_dummy("ode")
    model.solver_parameters = (dt=0.1, callback="test",)
    set_solver!(model, saveat=0.1, merge_kwargs=false)
    @test model.solver_parameters.saveat == 0.1  # added
    @test !(:dt in keys(model.solver_parameters))  # removed
    @test :callback in keys(model.solver_parameters)  # not removed

    # Throw an exception if the user attempts to pass a callback
    @test_throws OscillatorPopulationError set_solver!(model, callback="test")

end

@testset "get_parameter_index" begin

    model = create_dummy("ode")
    index = get_parameter_index(model, "c")
    @test index == 3

    # Nonexisting parameter
    @test_throws OscillatorPopulationError get_parameter_index(model, "d")

end

@testset "get_parameter" begin

    model = create_dummy("ode")
    value = get_parameter(model, "b")
    @test value == 2

    model = create_dummy("sde")
    value = get_parameter(model, "b")
    @test value == 2

    model = create_dummy("jump")
    value = get_parameter(model, "b")
    @test value == 2

    # Nonexisting parameter
    @test_throws OscillatorPopulationError get_parameter(model, "d")

end

@testset "set_parameter!" begin

    # Change one parameter
    model = create_dummy("ode")
    set_parameter!(model, "b", 123)
    @test model.problem.p[1] == 1  # not changed
    @test model.problem.p[2] == 123  # changed
    @test model.problem.p[3] == 0  # not changed

    # Change several parameters to the equal value
    model = create_dummy("sde")
    set_parameter!(model, "a=c", 777)
    @test model.problem.p[1] == 777  # changed
    @test model.problem.p[2] == 2  # not changed
    @test model.problem.p[3] == 777  # changed

    model = create_dummy("jump")
    set_parameter!(model, "c=b", -12)
    @test model.problem.prob.p[1] == 1  # changed
    @test model.problem.prob.p[2] == -12  # changed
    @test model.problem.prob.p[3] == -12  # changed

    # Try to change invalid parameter
    model = create_dummy("ode")
    @test_throws OscillatorPopulationError set_parameter!(model, "d", 0)

    # Change a list of parameters
    model = create_dummy("ode")
    set_parameter!(model, ["b"], [123])
    @test model.problem.p[1] == 1  # not changed
    @test model.problem.p[2] == 123  # changed
    @test model.problem.p[3] == 0  # not changed

    model = create_dummy("sde")
    set_parameter!(model, ["a", "c"], [111, 222])
    @test model.problem.p[1] == 111  # changed
    @test model.problem.p[2] == 2  # not changed
    @test model.problem.p[3] == 222  # changed

    model = create_dummy("jump")
    set_parameter!(model, ["a=c", "b"], [111, 222])
    @test model.problem.prob.p[1] == 111  # changed
    @test model.problem.prob.p[2] == 222  # changed
    @test model.problem.prob.p[3] == 111  # changed

    # Different array lengths for parameter names and values
    model = create_dummy("ode")
    @test_throws OscillatorPopulationError set_parameter!(model, ["a", "b"], [1, 2, 3])

    # Try to change invalid parameter in an array
    model = create_dummy("ode")
    @test_throws OscillatorPopulationError set_parameter!(model, ["a", "e"], [2, 3])

end

@testset "set_input!" begin

    model = create_dummy("ode")
    events = [2 5]
    set_input!(model, events, "b")

    # Check that the input field was updated
    @test model.input[1] == events
    @test model.input[2] == "b"

    # Check that a callback appeared in the solver parameters
    @test :callback in keys(model.solver_parameters)
    @test model.solver_parameters.callback isa CallbackSet

    # The original parameters are unchanged
    @test :dt in keys(model.solver_parameters)
    @test model.solver_parameters.dt == 1.0

    # Test that the callback can be rewritten
    callback_old = model.solver_parameters.callback
    @test callback_old == model.solver_parameters.callback  # not changed yet
    set_input!(model, events, "a")
    @test callback_old != model.solver_parameters.callback  # changed

    # For SDE a continuous callback is created
    model = create_dummy("sde")
    set_input!(model, events, "b")
    @test model.solver_parameters.callback isa CallbackSet

    # For a jump model a discrete callback is created
    model = create_dummy("jump")
    set_input!(model, events, "b")
    @test model.solver_parameters.callback isa DiscreteCallback

end

@testset "print_info" begin

    @test print_info isa Function

end

@testset "simulate_model" begin

    model = create_dummy("ode")
    x, t, sol = simulate_model(model)
    @test x == [1.0 0.0; 2.0 2.0; 3.0 4.0; 4.0 6.0; 5.0 8.0]
    @test t == [0.0, 1.0, 2.0, 3.0, 4.0]
    @test Matrix(sol[:, :]') == x
    @test sol.t == t

end