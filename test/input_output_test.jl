@testset "save_data" begin

    @test hasmethod(save_data, Tuple{DataFrame, Any}, (:force,))

end

@testset "load_data" begin

    @test hasmethod(load_data, Tuple{Any})

end

@testset "save_figure" begin

    @test hasmethod(save_figure, Tuple{Any, Any}, (:force,))

end
