using PowerSimulationsDynamicsSurrogates
using Documenter

DocMeta.setdocmeta!(
    PowerSimulationsDynamicsSurrogates,
    :DocTestSetup,
    :(using PowerSimulationsDynamicsSurrogates);
    recursive = true,
)

makedocs(;
    modules = [PowerSimulationsDynamicsSurrogates],
    authors = "Matt Bossart",
    repo = "https://github.com/m-bossart/PowerSimulationsDynamicsSurrogates.jl/blob/{commit}{path}#{line}",
    sitename = "PowerSimulationsDynamicsSurrogates.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://m-bossart.github.io/PowerSimulationsDynamicsSurrogates.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

# For running local web server of docs 
using LiveServer
serve(dir = "docs/build")

# Need to make repo public to use GitHub pages to deploy docs 
#= deploydocs(;
    repo = "github.com/m-bossart/PowerSimulationsDynamicsSurrogates.jl",
    devbranch = "main",
) =#
