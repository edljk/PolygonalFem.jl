#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-------------------------
"""
   IK, JK, SK, IM, JM, SM = assembKM_vemKM3D(pv, cellsb)

Assemble stiffness and mass terms of virtual element approach
"""
function assembKM_vemKM3D(pv, cellsb)
    # type of input (for overloaded differentiation calls)
    typein = typeof(pv[1])
    n_dofs, n_polys = size(pv, 1), 3 # method has 1 degree of freedom per vertex
    # stiffness coefficients and source coefficients   
    IK, JK = Int64[], Int64[]
    IM, JM = Int64[], Int64[]
    SK, SM = if typein == Float64 
        Float64[], Float64[]
    else
        Any[], Any[]
    end
    # impose an ordering on the linear polynomials
    linear_polynomials = [[0,0], [1,0], [0,1]] 
    H = zeros(typein, 3, 3)
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
        D = zeros(typein, n_sides, n_polys) 
        D[:, 1] .= 1
        B = zeros(typein, n_polys, n_sides) 
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
                D[vertex_id, poly_id] = dot(vert .- centroid[:], 
                                            poly_degree) / diameter
                B[poly_id, vertex_id] = 0.5 * dot(monomial_grad, vertex_normal)
            end
        end
        # routine to compute H
        fill!(H, zero(typein))
        H[1, 1] = area
        for j = 2:3
            for k = 2:3
                for v = 1:n_sides
                    vert = verts[v, :]
                    next = verts[mod_wrap(v - 1, n_sides), :]                
                    H[j, k] +=  _integrate2dtri(centroid, vert, next,
                      linear_polynomials[j], linear_polynomials[k], diameter, 6)
                end
            end
        end   
        # compute the local Ritz projector to polynomials
        projector = (B * D) \ B 
        stabilising_term = (eye(n_sides) .- D * projector)' * 
                           (eye(n_sides) .- D * projector)
        G = B * D 
        Gt = copy(G)
        Gt[1, :] .= 0
        local_stiffness = projector' * Gt * projector + stabilising_term
        C = H * (G \ B)
        ritz_projector = D * projector
        local_mass = C' * (H \ C) + area * (eye(n_sides) - ritz_projector)' *
                                           (eye(n_sides) - ritz_projector)
        # copy local to global
        append!(IK, repeat(vert_ids, length(vert_ids)))
        append!(JK, repeat(vert_ids', length(vert_ids))[:])
        append!(SK, local_stiffness[:])
        append!(IM, repeat(vert_ids, length(vert_ids)))
        append!(JM, repeat(vert_ids', length(vert_ids))[:])
        append!(SM, local_mass[:])
    end
    return IK, JK, SK, IM, JM, SM
end
function _integrate2dtri(P1, P2, P3, α, β, diameter, qrule)
    V = 0.
    if qrule == 1
        centroid = P1
        tricentroid = 1 / 3 * (P1 + P2 + P3)
        A = [1 P1...; 1 P2...; 1 P3...]
        area = 0.5 * abs(det(A))
        x, y = tricentroid[1], tricentroid[2]   
        ma = ((x - centroid[1]) / diameter) ^ α[1] *
             ((y - centroid[2]) / diameter) ^ α[2]
        mb = ((x - centroid[1]) / diameter) ^ β[1] *
             ((y - centroid[2]) / diameter) ^ β[2]
        V = area * ma * mb
    elseif qrule == 6
        centroid = P1
        A = [1 P1...; 1 P2...; 1 P3...]
        area = 0.5 * abs(det(A))  
        # 6th order triangular quadrature
        # obtain Quadrature points and weights    
        qw = [0.1116907948390, 0.1116907948390, 0.1116907948390, 
              0.0549758718276, 0.0549758718276, 0.0549758718276]
        qx = [0.108103018168, 0.445948490915, 0.445948490915, 
              0.816847572980, 0.091576213509, 0.091576213509]
        qy = [0.445948490915, 0.445948490915, 0.108103018168, 
              0.091576213509, 0.091576213509, 0.816847572980]
        V = 0.
        for q = 1:qrule
            xhat = (P2[1] - P1[1]) * qx[q] + (P3[1] - P1[1]) * qy[q] + P1[1]
            yhat = (P2[2] - P1[2]) * qx[q] + (P3[2] - P1[2]) * qy[q] + P1[2]
            ma = ((xhat - centroid[1]) / diameter) ^ α[1] *
                 ((yhat - centroid[2]) / diameter) ^ α[2]
            mb = ((xhat - centroid[1]) / diameter) ^ β[1] *
                 ((yhat - centroid[2]) / diameter) ^ β[2]
            V += qw[q] * 2 * area * ma * mb
        end
    end
    return V
end
#-------------------------------------------------------------------------------
"""
    u, p, t, meshboundary = vem_eigs(filename::String =  "squarepolmesh_coarse",
                                      nc::Int64 = 1_00; resolution::Int64 = 400)

N.B. 
* nc = 100, 1_000 or 10_000
* adapted code from the iVEM-1.0 Matlab Toolbox
"""
function vem_eigs(filename::String = "Lpolmesh", nc::Int64 = 1_00;
                  resolution::Int64 = 400, visible::Bool = false,
                  numeig::Int64 = 2)
    # computes the virtual element solution of the Poisson problem on 
    # the specified mesh
    # load the mesh  + pv       : vertices of the cells 
    #                + cells    : polygonal cells (⚠ orientation is crucial)
    #                + (pb, tb) : restricted delaunay mesh 
    mesh_filename = "$(@__DIR__)/../test/data/$(filename)_$(nc).jld2"
    println("read file $(mesh_filename)")
    JLD2.@load(mesh_filename, ps, Iv, pv, bbelem, elem, pb, tb) 
    cells, bcells = bbelem, bbelem
    # boundary points / rhs  
    boundary_condition, meshboundary = if occursin("square", mesh_filename)
        EPS = 1e-4
        dirichlet, union(findall(pv[:, 2] .> (1 - EPS)),
                         findall(pv[:, 2] .< EPS),
                         findall(pv[:, 1] .> (1 - EPS)),
                         findall(pv[:, 1] .< EPS))
    else
        EPS = 1e-4
        dirichlet, union(findall(pv[:, 2] .> (1 - EPS)),
                         findall(pv[:, 2] .< - 1 + EPS),
                         intersect(findall(abs.(pv[:, 2]) .< EPS),
                                   findall(pv[:, 1] .>= 0.)),
                         findall(pv[:, 1] .> (1 - EPS)),
                         findall(pv[:, 1] .< - 1 + EPS),
                         intersect(findall(abs.(pv[:, 1]) .< EPS),
                                   findall(pv[:, 2] .>= 0.)))
    end
    # call assemble function
    IK, JK, SK, IM, JM, SM = assembKM_vemKM(pv, cells)
    K, M, internal_dofs, u, Uu, λ = solve_eigs(IK, JK, SK, IM, JM, SM, pv,
                                               meshboundary, boundary_condition,
                                               numeig)
    @show round.(λ, digits = 6)
    # plot
    if resolution > 0
        Iv = visible ? Iv : 1:length(cells)
        plotsolution(u, pv, cells, bcells, Iv, ps,
                     mesh_filename = mesh_filename, resolution = resolution)
    end
    return nothing 
end
#-------------------------------------------------------------------------------
function solve_eigs(IK, JK, SK, IM, JM, SM, pv, meshboundary,
                    boundary_condition, numeig)
    n_dofs = size(pv, 1)   
    u = zeros(n_dofs) # degrees of freedom of the virtual element solution
    K = sparse(IK, JK, SK, n_dofs, n_dofs) 
    M = sparse(IM, JM, SM, n_dofs, n_dofs) 
    #@show maximum(abs.(K - K'))
    #@show maximum(abs.(M - M'))
    nb = length(meshboundary)
    boundary_vals = [boundary_condition(pv[meshboundary[k], :]) for k = 1:nb]
    # vertices which aren’t on the boundary
    internal_dofs =  setdiff(1:n_dofs, meshboundary) 
    # call eigs
    λ, Uu = eigs(K[internal_dofs, internal_dofs], 
                  M[internal_dofs, internal_dofs], nev = numeig, which = :LM, 
                  tol = 0., maxiter = 30_000, sigma = rand() / 1e4)  
    # solve
    u[internal_dofs] = real.(Uu[:, numeig])
    u[meshboundary] .= boundary_vals # set the boundary values
    return K, M, internal_dofs, u, Uu, real.(λ)
end
#-------------------------------------------------------------------------------
function dirichlet(points)
    return 0.
end