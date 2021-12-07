using SMLMFrameConnection
using Documenter

DocMeta.setdocmeta!(SMLMFrameConnection, :DocTestSetup, :(using SMLMFrameConnection); recursive=true)

makedocs(;
    modules=[SMLMFrameConnection],
    authors="klidke@unm.edu",
    repo="https://github.com/JuliaSMLM/SMLMFrameConnection.jl/blob/{commit}{path}#{line}",
    sitename="SMLMFrameConnection.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaSMLM.github.io/SMLMFrameConnection.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaSMLM/SMLMFrameConnection.jl",
)
