module PolygonalFem

using LinearAlgebra, StatsBase, SparseArrays, GroupSlices
using FileIO, JLD2, UnicodePlots

include("assembly.jl")
include("plot.jl")
eye(n) = Matrix(I, n, n)
"""
 adapted code from the article / matlab's code  "The virtual element method in 50 lines of matlab". See https://arxiv.org/pdf/1604.06021.pdf
"""
function vem(mesh_filename::String = 
                             "$(@__DIR__)/../test/data/squarepolmesh_100.jld2")
    # computes the virtual element solution of the Poisson problem on 
    # the specified mesh
    # load the mesh 
    JLD2.@load(mesh_filename, pv, cellsb, cellsbt) 
    nt = sum(map(x -> size(x, 1), cellsbt))
    t, ct = zeros(Int64, nt, 3), 1
    for k = 1:length(cellsbt)
        for l = 1:size(cellsbt[k], 1)
            t[ct, :] .= cellsbt[k][l, :]
            ct += 1
        end
    end
    meshboundary = unique(btri(t)[:])
    # identify boundary points
    n_dofs = size(pv, 1)
    n_polys = 3 # method has 1 degree of freedom per vertex
    K = spzeros(n_dofs, n_dofs) # stiffness matrix
    F = zeros(n_dofs, 1) # forcing vector
    u = zeros(n_dofs, 1) # degrees of freedom of the virtual element solution
    # impose an ordering on the linear polynomials
    linear_polynomials = [[0,0], [1,0], [0,1]] 
    mod_wrap(x, a) = mod(x - 1, a) + 1 # a utility function for wrapping around a vector
    areatot = 0.
    for el_id = 1:length(cellsb)
        vert_ids = cellsb[el_id]   # IDs of the vertices of this element
        verts = pv[vert_ids, :]    # coordinates of the vertices of this element
        n_sides = length(vert_ids) # start computing the geometric information
        area_components = verts[:, 1] .* [verts[2:end, 2]..., verts[1, 2]] .- 
                         [verts[2:end, 1]..., verts[2, 1]] .* verts[:, 2]
        area = 0.5 * abs(sum(area_components))
        areatot += area
        centroid = sum((verts .+ vcat(verts[2:end, :], verts[1, :]')) .* 
                        repeat(area_components, 1, 2), dims = 1) / 
                    (6 * area)
        diameter = 0 # compute the diameter by looking at every pair of vertices
        for i = 1:(n_sides-1)
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
            for poly_id = 2:n_polys #
                poly_degree = linear_polynomials[poly_id]
                # gradient of a linear polynomial is constant
                monomial_grad = poly_degree / diameter 
                D[vertex_id, poly_id] = dot(vert .- centroid[:], poly_degree) / diameter
                B[poly_id, vertex_id] = 0.5 * dot(monomial_grad, vertex_normal)
            end
        end
        # compute the local Ritz projector to polynomials
        projector = (B * D) \ B 
        stabilising_term = (eye(n_sides) - D * projector)' * (eye(n_sides) .- 
                            D * projector)
        G = B * D 
        G[1, :] .= 0
        local_stiffness = projector' * G * projector + stabilising_term
        # copy local to global
        K[vert_ids, vert_ids] .= K[vert_ids, vert_ids] .+ local_stiffness 
        F[vert_ids] .= F[vert_ids] .+ rhs(centroid) * area / n_sides
    end
    boundary_vals = boundary_condition(pv[meshboundary, :])
    # vertices which aren’t on the boundary
    internal_dofs =  setdiff(1:n_dofs, meshboundary) 
    # apply the boundary condition
    F -= K[:, meshboundary] * boundary_vals 
    # solve
    qrK = qr(K[internal_dofs, internal_dofs])
    u[internal_dofs] = qrK \ F[internal_dofs] 
    # K[internal_dofs, internal_dofs] \ F[internal_dofs] 
    u[meshboundary] .= boundary_vals # set the boundary values
    #plotunicode_solution(mesh, u)
    Main.closeall()
    f, L = Main.figure(1, fig3D = true)
    u ./= maximum(abs.(u))
    Main.plot_t(hcat(pv, u), t, scalars = u[:], representation = "surfacemesh", 
                colormap = "hsv", scene = L[1])
    Main.GraphicTools.finalizescene()
    return u
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