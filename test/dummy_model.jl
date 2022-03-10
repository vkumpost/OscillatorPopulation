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