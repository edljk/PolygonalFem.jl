module PolygonalFem

using LinearAlgebra, StatsBase, SparseArrays, GroupSlices
using FileIO, JLD2
using UnicodePlots, Makie, GeometryBasics, GLMakie, Colors, Interpolations

include("assembly.jl")
include("plot.jl")
eye(n) = Matrix(I, n, n)

"""
    u, p, t = vem(nc::Int64 = 1_00)

N.B. 
* nc = 100 or 1_000
* adapted code from the article / matlab's code  "The virtual element method in 50 lines of matlab". See https://arxiv.org/pdf/1604.06021.pdf
"""
function vem(nc::Int64 = 1_00)
    # computes the virtual element solution of the Poisson problem on 
    # the specified mesh
    # load the mesh  + pv       : vertices of the cells 
    #                + cellsb   : polygonal cells 
    #                + celssbt  : triangulated cells 
    #                + t        : all triangles 
    #                + (pb, tb) :  restricted delaunay mesh 
    mesh_filename = "$(@__DIR__)/../test/data/squarepolmesh_$(nc).jld2"
    println("read file $(mesh_filename)")
    JLD2.@load(mesh_filename, pv, cellsb, cellsbt, t, pb, tb) 
    # boundary points
    meshboundary = unique(btri(t)[:])
    # identify boundary points
    n_dofs, n_polys = size(pv, 1), 3 # method has 1 degree of freedom per vertex
    K = spzeros(n_dofs, n_dofs) # stiffness matrix
    F = zeros(n_dofs) # forcing vector
    u = zeros(n_dofs) # degrees of freedom of the virtual element solution
    # impose an ordering on the linear polynomials
    linear_polynomials = [[0,0], [1,0], [0,1]] 
    # a utility function for wrapping around a vector
    mod_wrap(x, a) = mod(x - 1, a) + 1 
    for el_id = 1:length(cellsb)
        #println(cellsb[el_id])
        vert_ids = cellsb[el_id]   # IDs of the vertices of this element
        verts = pv[vert_ids, :]    # coordinates of the vertices of this element
        n_sides = length(vert_ids) # start computing the geometric information
        area_components = verts[:, 1] .* [verts[2:end, 2]..., verts[1, 2]] .- 
                         [verts[2:end, 1]..., verts[1, 1]] .* verts[:, 2]
        area = 0.5 * abs(sum(area_components))
        centroid = sum((verts .+ vcat(verts[2:end, :], verts[1, :]')) .* 
                        repeat(area_components, 1, 2), dims = 1) / (6 * area)
        diameter = 0 # compute the diameter by looking at every pair of vertices
        for i = 1:(n_sides - 1)
            for j = (i + 1):n_sides
                diameter = max(diameter, norm(verts[i, :] .- verts[j, :]))
            end
        end
        D = zeros(n_sides, n_polys) 
        D[:, 1] .= 1
        B = zeros(n_polys, n_sides) 
        B[1, :] .= 1 / n_sides
        for vertex_id = 1:n_sides
            vert = verts[vertex_id, :] # this vertex and its neighbours
            prev = verts[mod_wrap(vertex_id - 1, n_sides), :]
            next = verts[mod_wrap(vertex_id + 1, n_sides), :]
            # average of normals on edges
            vertex_normal = [next[2] - prev[2], prev[1] - next[1]]
            # only need to loop over non-constant polynomials
            for poly_id = 2:n_polys 
                poly_degree = linear_polynomials[poly_id]
                # gradient of a linear polynomial is constant
                monomial_grad = poly_degree / diameter 
                D[vertex_id, poly_id] = dot(vert .- centroid[:], poly_degree) / diameter
                B[poly_id, vertex_id] = 0.5 * dot(monomial_grad, vertex_normal)
            end
        end
        # compute the local Ritz projector to polynomials
        projector = (B * D) \ B 
        stabilising_term = (eye(n_sides) .- D * projector)' * 
                           (eye(n_sides) .- D * projector)
        G = B * D 
        G[1, :] .= 0
        local_stiffness = projector' * G * projector + stabilising_term
        # copy local to global
        K[vert_ids, vert_ids] += local_stiffness 
        F[vert_ids] .+= rhs(centroid) * area / n_sides
    end
    boundary_vals = boundary_condition(pv[meshboundary, :])
    # vertices which aren’t on the boundary
    internal_dofs =  setdiff(1:n_dofs, meshboundary) 
    # apply the boundary condition
    F -= K[:, meshboundary] * boundary_vals 
    # solve
    u[internal_dofs] = K[internal_dofs, internal_dofs] \ F[internal_dofs] 
    u[meshboundary] .= boundary_vals # set the boundary values
    
    # plot
    plotsolution(u, pv, cellsb)
    return u, pv, t
end
function rhs(points)
    x, y = points[:, 1], points[:, 2]
    return 15 * sin.(π * x) .* sin.(π * y)
end
function boundary_condition(points)
    x, y = points[:, 1], points[:, 2]
    return (1 .- x) .* y .* sin.(π * x)
end

end # module