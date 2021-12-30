#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-------------------------
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