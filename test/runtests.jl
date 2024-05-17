using PlateauExperiment
using Test

@testset "PlateauExperiment.jl" begin
    
    @test τ_1(-1) == 0 && τ_1(1) == 1
    @test τ_2(0) == -1 && τ_2(1) == 1
    @test T_minus_one_one(-1; α = 4, β = 1) == -1 && T_minus_one_one(1; α = 4, β = 1) == -1
    @test T_minus_one_one(0; α = 4, β = 0.9) == 0.9
    
    @test all([T(x; α = 4, β = 0.9) == τ_1(T_minus_one_one(τ_2(x); α = 4, β = 0.9)) for x in 0:0.1:1])

end
