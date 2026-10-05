using PlateauExperiment
using Test

@testset "PlateauExperiment.jl" begin
    include("TestDynamic.jl") 
    include("TestOperators.jl") 
    include("TestAliasing.jl")
    include("TestLogDer.jl")
    include("TestExperiment.jl")  
end
