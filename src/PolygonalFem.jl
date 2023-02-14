module PolygonalFem

using LinearAlgebra, StatsBase, SparseArrays, GroupSlices
using NearestNeighbors, Arpack
using ReverseDiff, CoordinateTransformations, Interpolations
using ReverseDiff: GradientTape, GradientConfig, gradient, gradient!, compile
using FileIO, JLD2, GeometryBasics, Colors, PyCall, DelimitedFiles

include("assembly.jl")
include("assembly_elas.jl")
include("assembly_eigs.jl")
include("assembly_eigs3D.jl")
include("adjoint.jl")
include("utils.jl")
if ENV["HOSTNAME"] !== "ubuntuljk" || 
   :GraphicTools ∈ names(Main, all = true, imported = true)
    using  Makie, GLMakie
    include("plot.jl")
end
#-------------------------------------------------------------------------------
if :GraphicTools ∈ names(Main, all = true, imported = true)
    include("solveP1.jl")
    include("_genmeshes.jl")
end
eye(n) = Matrix(I, n, n)
#-------------------------------------------------------------------------------
"""
    u, p, t, meshboundary = vem(filename::String = "squarepolmesh_coarse",
                                nc::Int64 = 1_00;
                                resolution::Int64 = 400)

N.B. 
* nc = 100, 1_000 or 10_000
* adapted code from the article / matlab's code  "The virtual element method in 50 lines of matlab". See https://arxiv.org/pdf/1604.06021.pdf
"""
function vem(filename::String = "squarepolmesh_coarse", nc::Int64 = 1_00;
             resolution::Int64 = 400, visible::Bool = false)
    # computes the virtual element solution of the Poisson problem on 
    # the specified mesh
    # load the mesh  + pv       : vertices of the cells 
    #                + cells    : polygonal cells (⚠ orientation is crucial)
    #                + elem     : triangulated cells 
    #                + (pb, tb) : restricted delaunay mesh 
    mesh_filename = "$(@__DIR__)/../test/data/$(filename)_$(nc).jld2"
    println("read file $(mesh_filename)")
    JLD2.@load(mesh_filename, ps, Iv, pv, bbelem, elem, pb, tb) 
    cells, bcells = bbelem, bbelem
    # all triangles 
    t = vcat(elem...)
    # boundary points / rhs  
    rhs, boundary_condition = if occursin("square", mesh_filename)
        rhs_sqr, boundary_condition_sqr
    else
        rhs_L, boundary_condition_L
    end
    n_dofs, n_polys = size(pv, 1), 3 # method has 1 degree of freedom per vertex
    # identify boundary points
    meshboundary = unique(btri(t)[:])
    u = zeros(n_dofs) # degrees of freedom of the virtual element solution
    # call assemble function
    IK, JK, SK, IF, SF = assembKM_vemKsource(pv, cells, rhs,
                                             boundary_condition, meshboundary)
    K = sparse(IK, JK, SK, n_dofs, n_dofs) 
    F = Vector(sparsevec(IF, SF, n_dofs))
    boundary_vals = boundary_condition(pv[meshboundary, :])
    # vertices which aren’t on the boundary
    internal_dofs =  setdiff(1:n_dofs, meshboundary) 
    # solve
    u[internal_dofs] = K[internal_dofs, internal_dofs] \ F[internal_dofs] 
    u[meshboundary] .= boundary_vals # set the boundary values
    # plot
    if resolution > 0
        Iv = visible ? Iv : 1:length(cells)
        plotsolution(u, pv, cells, bcells, Iv, ps, 
                     mesh_filename = mesh_filename, resolution = resolution)
    end
    return u, pv, t, meshboundary
end
#-------------------------------------------------------------------------------
function rhs_sqr(points)
    x, y = points[:, 1], points[:, 2]
    return 15 * sin.(π * x) .* sin.(π * y)
end
#-------------------------------------------------------------------------------
function boundary_condition_sqr(points)
    x, y = points[:, 1], points[:, 2]
    return (1 .- x) .* y .* sin.(π * x)
    #return zeros(length(x))
end
#-------------------------------------------------------------------------------
function rhs_L(points)
    return zeros(size(points, 1))
end
#-------------------------------------------------------------------------------
function boundary_condition_L(points)
    x, y = points[:, 1], points[:, 2]
    pV = [[x[k], y[k]] for k = 1:length(x)]
    rθ = PolarFromCartesian().(pV)
    r, θ = [rθ[k].r for k = 1:length(x)],  [rθ[k].θ for k = 1:length(x)]
    return  r .^ (2 / 3.) .* sin.( 6 * (θ .- π / 2) / 3)
end
end # module