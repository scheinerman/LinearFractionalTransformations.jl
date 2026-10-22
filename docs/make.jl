using LinearFractionalTransformations
using Documenter

DocMeta.setdocmeta!(LinearFractionalTransformations, :DocTestSetup, :(using LinearFractionalTransformations); recursive=true)

makedocs(;
    modules=[LinearFractionalTransformations],
    authors="Ed Scheinerman and contributors",
    sitename="LinearFractionalTransformations.jl",
    checkdocs=:exports,
    format=Documenter.HTML(;
        canonical="https://scheinerman.github.io/LinearFractionalTransformations.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Construction" => "construction.md",
        "Evaluation" => "evaluation.md",
        "Operations" => "operations.md",
        "Infinity" => "infinity.md",
        "Projection" => "projection.md",
    ],
)

deploydocs(;
    repo="github.com/scheinerman/LinearFractionalTransformations.jl",
    devbranch="master",
)
