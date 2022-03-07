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