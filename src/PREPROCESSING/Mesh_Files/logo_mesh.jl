using Ferrite, FerriteGmsh, SparseArrays
using Downloads: download

function generate_mesh()
    logo_mesh = "logo.geo"
    asset_url = "https://raw.githubusercontent.com/Ferrite-FEM/Ferrite.jl/gh-pages/assets/"
    isfile(logo_mesh) || download(string(asset_url, logo_mesh), logo_mesh)

    grid = togrid(logo_mesh)

    addfacetset!(grid, "top",    x -> x[2] ≈ 1.0) # facets for which x[2] ≈ 1.0 for all nodes
    addfacetset!(grid, "left",   x -> abs(x[1]) < 1e-6)
    addfacetset!(grid, "bottom", x -> abs(x[2]) < 1e-6)

    return grid
end
