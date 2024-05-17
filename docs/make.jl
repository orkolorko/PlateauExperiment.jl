using PlateauExperiment
using Documenter

DocMeta.setdocmeta!(PlateauExperiment, :DocTestSetup, :(using PlateauExperiment); recursive=true)

makedocs(;
    modules=[PlateauExperiment],
    authors="Isaia Nisoli",
    sitename="PlateauExperiment.jl",
    format=Documenter.HTML(;
        canonical="https://orkolorko.github.io/PlateauExperiment.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/orkolorko/PlateauExperiment.jl",
    devbranch="main",
)
