using Test
using OpenPKPDCore

@testset "OpenPKPDCore skeleton" begin
    m = ModelSpec("test_model")
    s = SolverSpec("test_solver", 1e-6, 1e-9)
    result = simulate(m, s)

    @test result.metadata["deterministic"] == true
    @test result.model.name == "test_model"
    @test result.solver.reltol == 1e-6
end
