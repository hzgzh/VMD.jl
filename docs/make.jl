push!(LOAD_PATH,"../src/")
using Documenter, VMD

makedocs(
    modules = [VMD],
    sitename="VMD.jl",
    pages=Any[
        "Home" => "index.md",
        "Examples" => Dict(
        "example"  => "examples.md",
        "example1" => "example1.md",
        "example2" => "example2.md",
        "exmaple3" => "example3.md"),
        "Reference" => "api.md"
    ],
)



deploydocs(
    repo = "github.com/hzgzh/VMD.jl",
    branch = "gh-pages",
)