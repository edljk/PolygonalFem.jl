"""
    mes = meshmeasures(p, t)

compute areas of every element
"""
function meshareas(p, t)
    d12 = p[t[:, 2], :] - p[t[:, 1], :]
    d13 = p[t[:, 3], :] - p[t[:, 1], :]
    mes = (d12[:, 1] .* d13[:, 2] .- d12[:, 2] .* d13[:, 1]) / 2.
    return reshape(abs.(mes), length(mes))
end
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-------------------------
"""
   IK, JK, SK, IF, SF = assembKM_vemKsource(pv, cellsb, rhs, boundary_condition,
                                            meshboundary)

Assemble stiffness and source terms of virtual element approach
"""
function assembKM_vemKsource(pv, cellsb, rhs, boundary_condition,
                             meshboundary)
    # type of input (for overloaded differentiation calls)
    typein = typeof(pv[1])
    n_dofs, n_polys = size(pv, 1), 3 # method has 1 degree of freedom per vertex
    # stiffness coefficients and source coefficients   
    IK, JK, IF = Int64[], Int64[], Int64[]
    SK, SF = if typein == Float64 
        Float64[], Float64[]
    else
        Any[], Any[]
    end
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
        # K[vert_ids, vert_ids] += local_stiffness 
        # F[vert_ids] .+= rhs(centroid) * area / n_sides
        append!(IK, repeat(vert_ids, length(vert_ids)))
        append!(JK, repeat(vert_ids', length(vert_ids))[:])
        append!(SK, local_stiffness[:])
        append!(IF, vert_ids)
        append!(SF, repeat(rhs(centroid) * area / n_sides, length(vert_ids)))
    end
    # add boundary contribution in the source term 
    # old term F -= K[:, meshboundary] * boundary_vals 
    boundary_vals = boundary_condition(pv[meshboundary, :])
    Jinmesh = indexin(JK, meshboundary)
    Jboundary = findall(Jinmesh .!= nothing)
    for l ∈ Jboundary
        push!(IF, IK[l])
        push!(SF, - SK[l] * boundary_vals[Jinmesh[l]])
    end
    return IK, JK, SK, IF, SF
end
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-------------------------
""" 
   I, J = indKM_sparse(t)

non zero indices in stiffness and mass P1 sparse matrices
"""
function indKM_sparse(t)
    it1, it2, it3 = t[:, 1], t[:, 2], t[:, 3] 
    I = vcat(it1, it2, it3, it2, it3, it1, it1, it2, it3)
    J = vcat(it2, it3, it1, it1, it2, it3, it1, it2, it3)
    return I, J
end
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-------------------------
"""
     K, M = assembKM_P12D(p, t)

assemble the P1 stifness and mass matrices.
"""
function assembKM_P12D(p, t)
    np, dim, nt = size(p, 1), size(p, 2), size(t, 1)
    # compute ar add define arK and arM (not used ?)
    ar = meshareas(p, t)
    arK, arM = ar, ar
    it1, it2, it3 = t[:, 1], t[:, 2], t[:, 3] 
    g1x, g1y, g2x, g2y, g3x, g3y = pdetrg(p, t, ar)
    c3 = (g1x .* g2x .+ g1y .* g2y) .* arK
    c1 = (g2x .* g3x .+ g2y .* g3y) .* arK
    c2 = (g3x .* g1x .+ g3y .* g1y) .* arK
    SK = vcat(c3, c1, c2, c3, c1, c2, -c2 - c3, -c3 - c1, -c1 - c2)
    # mass matrix
    aod = arM / 12. # off diagonal element
    ad = 2 * aod    # diagonal element
    SM = vcat(aod, aod, aod, aod, aod, aod, ad, ad, ad)
    # matrices 
    I, J = indKM_sparse(t)
    K = sparse(I, J, SK, np, np)
    M = sparse(I, J, SM, np, np)
    return K, M
end 
#-------------------------------------------------------------------------------
"""
   g1x, g1y, g2x, g2y, g3x, g3y = pdetrg(p, t, ar)

returns the  the gradient components of the triangle base functions
"""
function pdetrg(p, t, ar)
    # corner point indices
    a1, a2, a3 =  t[:, 1], t[:, 2], t[:, 3]
    # triangle sides
    r23x =  p[a3, 1] .- p[a2, 1]
    r23y =  p[a3, 2] .- p[a2, 2]
    r31x =  p[a1, 1] .- p[a3, 1]
    r31y =  p[a1, 2] .- p[a3, 2]
    r12x =  p[a2, 1] .- p[a1, 1]
    r12y =  p[a2, 2] .- p[a1, 2]    
    g1x =  - 0.5 .* r23y ./ ar
    g1y =    0.5 .* r23x ./ ar
    g2x =  - 0.5 .* r31y ./ ar
    g2y =    0.5 .* r31x ./ ar
    g3x =  - 0.5 .* r12y ./ ar
    g3y =    0.5 .* r12x ./ ar
    return g1x, g1y, g2x, g2y, g3x, g3y
end
#-------------------------------------------------------------------------------
"""
    be = btri(T::Array{Int64, 2})

finds boundary edges
"""
function btri(T::Array{Int64, 2})
    # form all faces, non - duplicates are surface triangles
    faces = vcat(T[:, [1, 2]], T[:, [1, 3]], T[:, [2, 3]])
    faces = sort(faces, dims = 2)
    _, ix, jx = uniquerows_int(faces)
    ee = collect(1:(maximum(jx) + 1)) .- 1e-8
    result = fit(Histogram, jx, (1:(maximum(jx) + 1)) .- 1e-8, 
                 closed = :right)
    vec = result.weights
    qx = findall(x -> x != 0, vec .- 2) 
    return faces[ix[qx], :]
end
#-------------------------------------------------------------------------------
function uniquerows_int(x::Matrix{Int64})
    Vx = [x[k, :] for k = 1:size(x, 1)]
    jj = groupslices(Vx) 
    ix, jx = firstinds(jj), similar(jj)
    for k = 1:length(ix)
        jx[jj .== ix[k]] .= k 
    end
    Vy = Vx[ix]
    y = [Vy[k][l] for k = 1:length(Vy), l = 1:length(Vy[1])]
    return y, ix, jx
end