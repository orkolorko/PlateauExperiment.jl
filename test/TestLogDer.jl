α = 3.5
β = 1
ln = PlateauExperiment.FourierLogDer(α, β, 128)

vals = [-0.5540898;
        1.47372468;
        -0.5642645]

@test vals[1] * ln[1] >= 0
@test vals[2] * ln[2] >= 0
@test vals[3] * ln[3] >= 0
@test vals[2] * ln[end] >= 0
@test vals[3] * ln[end - 1] >= 0
