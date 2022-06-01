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
function costnormU(pvin, cellsb, t, rhs, boundary_condition, meshboundary)
    # copy coordinates
    pv = reshape(pvin, length(pvin) ÷ 2, 2)
    n_dofs, n_polys = size(pv, 1), 3 # method has 1 degree of freedom per vertex
    u = zeros(n_dofs) # degrees of freedom of the virtual element solution
    # call assemble function
    IK, JK, SK, IF, SF = assembKM_vemKsource(pv, cellsb, rhs,
                                             boundary_condition, meshboundary)
    K = sparse(IK, JK, SK, n_dofs, n_dofs) 
    F = Vector(sparsevec(IF, SF, n_dofs))
    boundary_vals = boundary_condition(pv[meshboundary, :])
    # vertices which aren’t on the boundary
    internal_dofs =  setdiff(1:n_dofs, meshboundary) 
    # solve
    u[internal_dofs] = K[internal_dofs, internal_dofs] \ F[internal_dofs] 
    u[meshboundary] .= boundary_vals # set the boundary values
    # L² cost 
    val = sum(u[k] ^ 2 for k = 1:n_dofs) / n_dofs
    return val
end
#-------------------------------------------------------------------------------
function ∇costnormU(pvin, cellsb, t, rhs, boundary_condition, meshboundary)
    # copy coordinates
    pv = reshape(pvin, length(pvin) ÷ 2, 2)
    n_dofs, n_polys = size(pv, 1), 3 # method has 1 degree of freedom per vertex
    u = zeros(n_dofs) # degrees of freedom of the virtual element solution
    # call assemble function
    IK, JK, SK, IF, SF = assembKM_vemKsource(pv, cellsb, rhs,
                                             boundary_condition, meshboundary)
    K = sparse(IK, JK, SK, n_dofs, n_dofs) 
    F = Vector(sparsevec(IF, SF, n_dofs))
    boundary_vals = boundary_condition(pv[meshboundary, :])
    # vertices which aren’t on the boundary
    internal_dofs =  setdiff(1:n_dofs, meshboundary) 
    # solve
    u[internal_dofs] = K[internal_dofs, internal_dofs] \ F[internal_dofs] 
    u[meshboundary] .= boundary_vals # set the boundary values
    # cost gradient and adjoint function
    JFU = - 2 * u / n_dofs
    λ = K' \ JFU
    λ = zeros(length(u))
    λ[internal_dofs] = K'[internal_dofs, internal_dofs] \ JFU[internal_dofs]  
    # compute full gradient by automatic differentiation
    fprod(x) = prodgrad(x, cellsb, t, rhs, boundary_condition, meshboundary, 
                        u, λ)
    ∇prodgrad = ReverseDiff.gradient(fprod, pv)
    #= optimization tests
    ∇prodgrad2 = similar(∇prodgrad)
    f_tape = GradientTape(fprod, randn(size(pv)))
    compiled_f_tape = compile(f_tape)
    @time gradient!(∇prodgrad2, compiled_f_tape, pv)
    println("diff grads = $(norm(∇prodgrad2 - ∇prodgrad2))")
    =#
    # gradient contribution of boundary points
    fboundary(x) = sum(boundary_condition(x[meshboundary, :]) .^ 2) / n_dofs
    ∇Dirichlet = ReverseDiff.gradient(fboundary, pv)
    return ∇prodgrad + ∇Dirichlet
end
#-------------------------------------------------------------------------------
function prodgrad(pv, cellsb, t, rhs, boundary_condition, meshboundary, u, λ)
    # assemble matrices
    IK, JK, SK, IF, SF = assembKM_vemKsource(pv, cellsb, rhs,
                                             boundary_condition, meshboundary)
    val = 0.
    # consider only internal points
    Jinmesh = indexin(JK, meshboundary)
    Jboundary = findall(Jinmesh .!= nothing)
    SK[Jboundary] .= 0.
    for (i, j, s) ∈ zip(IK, JK, SK)
        val += λ[i] * s * u[j]
    end
    for (i, s) ∈ zip(IF, SF) 
        val -= λ[i] * s 
    end
    return val
end
#-------------------------------------------------------------------------------
function costnormU_elas2D(pvin, cellsb, t, rhs, boundary_condition,
                          meshboundary)
    # copy coordinates
    pv = reshape(pvin, length(pvin) ÷ 2, 2)
    # call assemble function
    IK, JK, SK, IF, SF = assembKM_vemKsource_elas2D(pv, cellsb, rhs,
                                                     boundary_condition, 
                                                     meshboundary)
    # solve 
    K, F, internal_dofs, u = solve_elas2D(IK, JK, SK, IF, SF, pv, meshboundary, 
                                          boundary_condition)
    # L² cost 
    n_dofs = length(u)
    val = sum(u[k] ^ 2 for k = 1:n_dofs) / n_dofs
    return val
end
#-------------------------------------------------------------------------------
function ∇costnormU_elas2D(pvin, cellsb, t, rhs, boundary_condition, 
                           meshboundary)
    # copy coordinates
    pv = reshape(pvin, length(pvin) ÷ 2, 2)
    # call assemble function
    IK, JK, SK, IF, SF = assembKM_vemKsource_elas2D(pv, cellsb, rhs,
                                             boundary_condition, meshboundary)
    # solve 
    K, F, internal_dofs, u = solve_elas2D(IK, JK, SK, IF, SF, pv, meshboundary, 
                                          boundary_condition)
    n_dofs = length(u)
    # cost gradient and adjoint function
    JFU = - 2 * u / n_dofs
    λ = K' \ JFU
    λ = zeros(n_dofs)
    λ[internal_dofs] = K'[internal_dofs, internal_dofs] \ JFU[internal_dofs]  
    # compute full gradient by automatic differentiation
    fprod(x) = prodgrad_elas2D(x, cellsb, t, rhs, boundary_condition,
                               meshboundary, u, λ)
    ∇prodgrad = ReverseDiff.gradient(fprod, pv)
    # gradient contribution of boundary points
    fboundary(x) = sum(boundary_condition(x[meshboundary, :]) .^ 2) / n_dofs
    ∇Dirichlet = ReverseDiff.gradient(fboundary, pv)
    return ∇prodgrad #+ ∇Dirichlet
end
#-------------------------------------------------------------------------------
function prodgrad_elas2D(pv, cellsb, t, rhs, boundary_condition, 
                         meshboundary, u, λ)
    # assemble matrices
    IK, JK, SK, IF, SF = assembKM_vemKsource_elas2D(pv, cellsb, rhs,
                                                    boundary_condition, 
                                                    meshboundary)
    val = 0.
    # boundary dofs
    #=
    Vmeshboundary = vcat(2 * meshboundary .- 1, 2 * meshboundary)
    Jinmesh = indexin(JK, Vmeshboundary)
    Jboundary = findall(Jinmesh .!= nothing)
    SK[Jboundary] .= 0.
    =#
    for (i, j, s) ∈ zip(IK, JK, SK)
        val += λ[i] * s * u[j]
    end
    for (i, s) ∈ zip(IF, SF) 
        val -= λ[i] * s 
    end
    return val
end
#-------------------------------------------------------------------------------
function costeigs2D(pvin, cellsb, boundary_condition, meshboundary, numeig)
    # copy coordinates
    pv = reshape(pvin, length(pvin) ÷ 2, 2)
    # call assemble function
    IK, JK, SK, IM, JM, SM = assembKM_vemKM(pv, cellsb)
    # compute eigenvalues
    K, M, internal_dofs, u, Uu, λ = solve_eigs(IK, JK, SK, IM, JM, SM, pv,
                                               meshboundary, boundary_condition,
                                               numeig)
    return λ[numeig], u
end
#-------------------------------------------------------------------------------
function ∇costeigs2D(pvin, cellsb, boundary_condition, meshboundary, numeig)
    # copy coordinates
    pv = reshape(pvin, length(pvin) ÷ 2, 2)
    n_dofs = size(pv, 1)   
    # call assemble function
    # call assemble function
    IK, JK, SK, IM, JM, SM = assembKM_vemKM(pv, cellsb)
    M = sparse(IM, JM, SM, n_dofs, n_dofs) 
    # compute eigenvalues
    K, M, internal_dofs, u, Uu, λ = solve_eigs(IK, JK, SK, IM, JM, SM, pv,
                                       meshboundary, boundary_condition, numeig)
    n_dofs = length(u)
    # compute full gradient by automatic differentiation
    fprod(x) = prodgrad_eigs2D(x, cellsb, boundary_condition, meshboundary,
                               u, λ[numeig])
    ∇prodgrad = ReverseDiff.gradient(fprod, pv) / (u' * M * u)
    # gradient contribution of boundary points
    #fboundary(x) = sum(boundary_condition(x[meshboundary, :]) .^ 2) / n_dofs
    #∇Dirichlet = ReverseDiff.gradient(fboundary, pv)
    return ∇prodgrad #+ ∇Dirichlet
end
#-------------------------------------------------------------------------------
function prodgrad_eigs2D(pv, cellsb, boundary_condition,  meshboundary, u, λ)
    # assemble matrices
    IK, JK, SK, IM, JM, SM = assembKM_vemKM(pv, cellsb)
    valK, valM = 0., 0.
    # boundary dofs
    #=
    Vmeshboundary = vcat(2 * meshboundary .- 1, 2 * meshboundary)
    Jinmesh = indexin(JK, Vmeshboundary)
    Jboundary = findall(Jinmesh .!= nothing)
    SK[Jboundary] .= 0.
    =#
    for (i, j, s) ∈ zip(IK, JK, SK)
        valK += u[i] * s * u[j] 
    end
    for (i, j, s) ∈ zip(IM, JM, SM)
        valM += u[i] * s * u[j] 
    end
    return valK - λ * valM
end