@testset "save_data, load_data, save_figure" begin

    @test save_data isa Function
    @test load_data isa Function
    @test save_figure isa Function

end
