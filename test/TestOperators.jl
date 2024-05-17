@testset "Operators" begin
    D = NoiseInterval(1, 2)

    @test D[1, 1] == 1
    @test exp(-π^2 / 2) ∈ D[2, 2]
    @test exp(-π^2 / 2) ∈ D[5, 5]

    D = NoiseBall(1, 2)
    @test 1 ∈ D[1, 1]
    @test exp(-π^2 / 2) ∈ D[2, 2]
    @test exp(-π^2 / 2) ∈ D[5, 5]

    v = [1; im; im; 0; 0]
    @test symmetrize_density(v) == [1; im; im; -im; -im]

end