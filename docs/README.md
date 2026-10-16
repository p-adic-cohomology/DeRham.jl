# Building the documentation

From the repository root, instantiate the documentation environment and build
the site:

```sh
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The generated site is written to `docs/build/`. The documentation workflow
builds and deploys the site on pushes to `main` and version tags, and builds it
for pull requests. Documenter publishes previews for pull requests from this
repository; previews are unavailable for fork pull requests because their
workflow tokens cannot write to the deployment branch. Closed pull requests
have their preview removed from `gh-pages` by the preview-cleanup workflow.

For a live preview, add LiveServer to the documentation environment if needed,
then run `using LiveServer; servedocs()` from the repository root.
