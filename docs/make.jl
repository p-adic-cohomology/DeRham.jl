using DeRham
using Documenter

makedocs(
    sitename = "DeRham.jl",
    modules = [DeRham],
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting_started.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        edit_link = "main",
    ),
)

deploydocs(
    repo = "github.com/UCSD-computational-number-theory/DeRham.jl.git",
    devbranch = "main",
    push_preview = true,
)
