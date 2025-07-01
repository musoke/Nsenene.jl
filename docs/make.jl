using Nsenene
using Documenter

DocMeta.setdocmeta!(Nsenene, :DocTestSetup, :(using Nsenene); recursive=true)

makedocs(;
    modules=[Nsenene],
    authors="Nathan Musoke <nathan.musoke@gmail.com> and contributors",
    sitename="Nsenene.jl",
    format=Documenter.HTML(;
        canonical="https://musoke.github.io/Nsenene.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/musoke/Nsenene.jl",
    devbranch="main",
)
