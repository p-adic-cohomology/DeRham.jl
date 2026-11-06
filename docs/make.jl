using DeRham
using Documenter

function readme_examples()
    readme = read(joinpath(@__DIR__, "..", "README.md"), String)
    opening = "<!-- examples -->"
    closing = "<!-- /examples -->"
    sections = split(readme, opening; limit = 2)
    length(sections) == 2 || error("README is missing the examples opening marker")
    examples = split(sections[2], closing; limit = 2)
    length(examples) == 2 || error("README is missing the examples closing marker")
    return strip(examples[1])
end

template = read(joinpath(@__DIR__, "templates", "getting_started.md"), String)
generated = joinpath(@__DIR__, "src", "getting_started.md")
write(generated, string(rstrip(template), "\n\n## Options\n\n", readme_examples(), "\n"))

DocMeta.setdocmeta!(DeRham, :DocTestSetup, :(using DeRham, Oscar); recursive = true)

makedocs(
    sitename = "DeRham.jl",
    modules = [DeRham],
    checkdocs = :none,
    doctest = true,
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting_started.md",
        "Invariants" => [
            "Zeta function" => "invariants/zeta.md",
            "Frobenius matrix" => "invariants/frobenius_matrix.md",
            "Newton & Hodge polygons" => "invariants/polygons.md",
            "Hasse–Witt matrix & a-number" => "invariants/hasse_witt.md",
            "Picard rank" => "invariants/picard_rank.md",
            "Artin–Mazur height" => "invariants/artin_mazur_height.md",
            "Artin invariant" => "invariants/artin_invariant.md",
            "Jacobians" => "invariants/jacobians.md",
        ],
        "Utilities" => "utilities.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        edit_link = "main",
    ),
)

deploydocs(
    repo = "github.com/p-adic-cohomology/DeRham.jl.git",
    devbranch = "main",
    push_preview = true,
)
