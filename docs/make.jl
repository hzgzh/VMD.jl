push!(LOAD_PATH,"../src/")
using Documenter, VMD

makedocs(
    modules = [VMD],
    sitename="VMD.jl",
    pages=Any[
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Reference" => "api.md"
    ],
)



deploydocs(
    repo = "github.com/hzgzh/VMD.jl",
    branch = "gh-pages",
)