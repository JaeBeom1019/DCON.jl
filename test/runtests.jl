using Test
using DCON

@testset "DCON loads and runs basic" begin
    ini = DCON.Input()
    gfile_path = joinpath(@__DIR__, "..", "g000000.00000")  # 상대경로 안전하게
    
    dcon_result, δW = DCON.run_dcon(ini, gfile_path)
    
    @test isa(dcon_result, Dict)
    @test !isempty(δW)
end