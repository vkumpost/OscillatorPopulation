@testset "load_model" begin
    
    @test_throws OscillatorPopulationError load_model("nonexisting model")

end

@testset "amplitude-phase model" begin
    
    # ode model
    model = load_model("amplitude-phase", "ode")
    @test model.variable_names == ["x", "y"]
    @test model.parameter_names == ["λ", "A", "T", "I"]
    @test model.problem isa ODEProblem
    @test model.problem.u0 == [0.1, 0.1]
    @test model.problem.tspan == (0.0, 10.0)
    @test model.problem.p == [0.5, 2.0, 1.0, 0.0]
    @test model.solver_algorithm == DP5()
    @test model.solver_parameters == (saveat=0.01, reltol=1e-9, abstol=1e-9,)
    @test isempty(model.input[1])
    @test model.input[1] isa Matrix
    @test isempty(model.input[2])
    @test model.input[2] isa String
    @test model.output isa Function

    # sde model
    model = load_model("amplitude-phase", "sde")
    @test model.variable_names == ["x", "y"]
    @test model.parameter_names == ["λ", "A", "T", "I", "σ"]
    @test model.problem isa SDEProblem
    @test model.problem.u0 == [0.1, 0.1]
    @test model.problem.tspan == (0.0, 10.0)
    @test model.problem.p == [0.5, 2.0, 1.0, 0.0, 0.5]
    @test model.solver_algorithm == SOSRI()
    @test model.solver_parameters == (saveat=0.01,)
    @test isempty(model.input[1])
    @test model.input[1] isa Matrix
    @test isempty(model.input[2])
    @test model.input[2] isa String
    @test model.output isa Function

    # jump model
    @test_throws OscillatorPopulationError load_model("amplitude-phase", "jump")

end

@testset "goodwin model" begin

    # ode model
    model = load_model("goodwin", "ode")
    @test model.variable_names == ["x", "y", "z"]
    @test model.parameter_names == ["K", "n", "a1", "a2", "a3", "d1", "d2", "d3"]
    @test model.problem isa ODEProblem
    @test model.problem.u0 == [0.1, 0.1, 0.1]
    @test model.problem.tspan == (0.0, 100.0)
    @test model.problem.p == [1.0, 10.0, 5.0, 5.0, 5.0, 0.5, 0.5, 0.5]
    @test model.solver_algorithm == DP5()
    @test model.solver_parameters == (saveat=0.01, reltol=1e-9, abstol=1e-9,)
    @test isempty(model.input[1])
    @test model.input[1] isa Matrix
    @test isempty(model.input[2])
    @test model.input[2] isa String
    @test model.output isa Function

    # sde model
    @test_throws OscillatorPopulationError load_model("goodwin", "sde")

    # jump model
    @test_throws OscillatorPopulationError load_model("goodwin", "jump")
    
end

@testset "goodwin-general model" begin

    # ode model
    model = load_model("goodwin-general", "ode")
    @test model.variable_names == ["x1", "x2", "x3"]
    @test model.parameter_names == ["K", "n"]
    @test model.problem isa ODEProblem
    @test model.problem.u0 == [0.1, 0.1, 0.1]
    @test model.problem.tspan == (0.0, 10.0)
    @test model.problem.p == [0.1, 12.0]
    @test model.solver_algorithm == DP5()
    @test model.solver_parameters == (saveat=0.01, reltol=1e-9, abstol=1e-9,)
    @test isempty(model.input[1])
    @test model.input[1] isa Matrix
    @test isempty(model.input[2])
    @test model.input[2] isa String
    @test model.output isa Function

    model = load_model("goodwin-general", "ode", n_equations=2)
    @test model.variable_names == ["x1", "x2"]
    @test model.parameter_names == ["K", "n"]
    @test model.problem isa ODEProblem
    @test model.problem.u0 == [0.1, 0.1]
    @test model.problem.tspan == (0.0, 10.0)
    @test model.problem.p == [0.1, 12.0]
    @test model.solver_algorithm == DP5()
    @test model.solver_parameters == (saveat=0.01, reltol=1e-9, abstol=1e-9,)
    @test isempty(model.input[1])
    @test model.input[1] isa Matrix
    @test isempty(model.input[2])
    @test model.input[2] isa String
    @test model.output isa Function

    # sde model
    @test_throws OscillatorPopulationError load_model("goodwin-general", "sde")

    # jump model
    @test_throws OscillatorPopulationError load_model("goodwin-general", "jump")
    
end

@testset "van-der-pol model" begin
    
    # ode model
    model = load_model("van-der-pol", "ode")
    @test model.variable_names == ["x", "y"]
    @test model.parameter_names == ["B", "d", "I"]
    @test model.problem isa ODEProblem
    @test model.problem.u0 == [0.1, 0.1]
    @test model.problem.tspan == (0.0, 100.0)
    @test model.problem.p == [10.0, 2.0, 0.0]
    @test model.solver_algorithm == DP5()
    @test model.solver_parameters == (saveat=0.01, reltol=1e-9, abstol=1e-9,)
    @test isempty(model.input[1])
    @test model.input[1] isa Matrix
    @test isempty(model.input[2])
    @test model.input[2] isa String
    @test model.output isa Function

    # sde model
    model = load_model("van-der-pol", "sde")
    @test model.variable_names == ["x", "y"]
    @test model.parameter_names == ["B", "d", "I", "σ"]
    @test model.problem isa SDEProblem
    @test model.problem.u0 == [0.1, 0.1]
    @test model.problem.tspan == (0.0, 100.0)
    @test model.problem.p == [10.0, 2.0, 0.0, 0.1]
    @test model.solver_algorithm == SOSRI()
    @test model.solver_parameters == (saveat=0.01,)
    @test isempty(model.input[1])
    @test model.input[1] isa Matrix
    @test isempty(model.input[2])
    @test model.input[2] isa String
    @test model.output isa Function

    # jump model
    @test_throws OscillatorPopulationError load_model("van-der-pol", "jump")

end
