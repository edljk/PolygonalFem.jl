"""
    fill_vertices!(pvector, parray)
"""
function fill_vertices!(pvector, parray)
    dim = size(parray[1], 1)
    # check dimensions 
    if sum(length.(parray)) != length(pvector)
        error("pb of dimensions in fill_vertices $(sum(length.(parray))) != $(length(pvector)) ")
    end
    c = 1
    for k = 1:length(parray)
        for l = 1:dim # fill by column 
            for m = 1:size(parray[k], 1)
                parray[k][m, l] = pvector[c]
                c += 1
            end
        end
    end
    return nothing 
end
#-------------------------------------------------------------------------------
function costnormU(pvin, cellsb, t, rhs, boundary_condition)
    # copy coordinates
    pv = reshape(pvin, length(pvin) ÷ 2, 2)
    n_dofs, n_polys = size(pv, 1), 3 # method has 1 degree of freedom per vertex
    # identify boundary points
    meshboundary = unique(btri(t)[:])
    u = zeros(n_dofs) # degrees of freedom of the virtual element solution
    # call assemble function
    IK, JK, SK, IF, SF = assembKM_vemKsource(pv, cellsb, rhs)
    K = sparse(IK, JK, SK, n_dofs, n_dofs) 
    F = Vector(sparsevec(IF, SF, n_dofs))
    boundary_vals = boundary_condition(pv[meshboundary, :])
    # vertices which aren’t on the boundary
    internal_dofs =  setdiff(1:n_dofs, meshboundary) 
    # apply the boundary condition
    F -= K[:, meshboundary] * boundary_vals 
    # solve
    u[internal_dofs] = K[internal_dofs, internal_dofs] \ F[internal_dofs] 
    u[meshboundary] .= boundary_vals # set the boundary values
    # L² cost 
    val = sum(u[k] ^ 2 for k = 1:n_dofs) / n_dofs
    return val
end
#-------------------------------------------------------------------------------
function ∇costnormU(pvin, cellsb, t, rhs, boundary_condition)
    # copy coordinates
    pv = reshape(pvin, length(pvin) ÷ 2, 2)
    n_dofs, n_polys = size(pv, 1), 3 # method has 1 degree of freedom per vertex
    # identify boundary points
    meshboundary = unique(btri(t)[:])
    u = zeros(n_dofs) # degrees of freedom of the virtual element solution
    # call assemble function
    IK, JK, SK, IF, SF = assembKM_vemKsource(pv, cellsb, rhs)
    K = sparse(IK, JK, SK, n_dofs, n_dofs) 
    F = Vector(sparsevec(IF, SF, n_dofs))
    boundary_vals = boundary_condition(pv[meshboundary, :])
    # vertices which aren’t on the boundary
    internal_dofs =  setdiff(1:n_dofs, meshboundary) 
    # apply the boundary condition
    F -= K[:, meshboundary] * boundary_vals 
    # solve
    u[internal_dofs] = K[internal_dofs, internal_dofs] \ F[internal_dofs] 
    u[meshboundary] .= boundary_vals # set the boundary values
      # cost gradient and adjoint function
    JFU = - 2 * u / n_dofs
    λ = K' \ JFU
    # compute full gradient by automatic differentiation
    ∇prodgrad = ReverseDiff.gradient(x -> prodgrad(x, cellsb, t, rhs, u, λ), pv)
    return ∇prodgrad
end
#-------------------------------------------------------------------------------
function prodgrad(pv, cellsb, t, rhs, u, λ)
    # assemble matrices
    IK, JK, SK, IF, SF = assembKM_vemKsource(pv, cellsb, rhs)
    val = 0.
    for (i, j, s) ∈ zip(IK, JK, SK)
        val += λ[i] * s * u[j]
    end
    for (i, s) ∈ zip(IF, SF)
        val -= λ[i] * s 
    end
    return val
end