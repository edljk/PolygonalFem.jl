include("assembly_eigs3D_utils.jl")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-------------------------
"""
    assembKM_vemKM3D(pv, elem, matbbelem, c)  → IK, JK, SK, IM, JM, SM 

Assemble stiffness and mass terms of virtual element approach.
The code is a julia translation of [Yu yue](https://github.com/Terenceyuyue/mVEM/tree/master/vem3)

"""
function assembKM_vemKM3D(pv::Matrix{T}, elem, matbbelem, c) where T
    # type of input (for overloaded differentiation calls)
    nps, typein = length(elem), typeof(pv[1])
    # compute centroids / diameters
    ps = zeros(nps)
    _, centroids, _, diameters = centers_momentum_geogram(ps, pv, elem)
    # define face and elem2face
    nbf_by_elem = length.(matbbelem)
    allfaces, elem2face = Vector{Vector{Int64}}(), Vector{Vector{Int64}}()
    for k = 1:nps 
        append!(allfaces, deepcopy(matbbelem[k]))
    end
    allfaces_keep = deepcopy(allfaces)
    sort!.(allfaces)
    ordered_faces = unique(allfaces)
    Ioldorder = [findfirst(isequal(ff), allfaces) for ff ∈ ordered_faces]
    faces = allfaces_keep[Ioldorder, :]
    If = [findfirst(isequal(allfaces[l]), ordered_faces) 
                                                    for l = 1:length(allfaces)]
    ck = 1
    for k = 1:nps 
        push!(elem2face, If[ck:(ck + nbf_by_elem[k] - 1)])
        ck += nbf_by_elem[k]
    end
    N, NT, NF = size(pv, 1), size(elem, 1), length(faces)
    # derive elliptic projections of all faces
    faceProj = Matrix{T}[]
    for (k, f) ∈ enumerate(faces)
        # pifs
        Pifs = faceEllipticProjection(pv[f, :])
        # sort columns
        idx = sortperm(f)
        push!(faceProj, copy(Pifs[:, idx]))
    end
    # compute and assemble the linear system
    Ph = Matrix{T}[] # matrix for error evaluation
    elem2dof = Vector{Int64}[]
    elemLen = [length(unique(union(matbbelem[k]...))) for k = 1:nps]
    nnz = sum(elemLen .^ 2)
    ii, jj = zeros(Int64, nnz), zeros(Int64, nnz)
    SK, SM = zeros(T, nnz), zeros(T, nnz)
    nnz = sum(elemLen)
    elemb = zeros(Int64, nnz)
    ia, ib = 0, 0
    for iel = 1:NT
        # ------- element information --------
        # faces
        elemf, indexFace = matbbelem[iel], elem2face[iel]
        # global index of vertices and local index of elemf
        Tri, indexDof, ~, elemfLocal = faceTriangulation(elemf)
        # centroid and diameter
        Nv = length(indexDof)
        Ndof = Nv
        V = @view pv[indexDof, :]
        ndb = 20
        #println("^" ^ ndb, " using approx centroid for debugging! ", "^" ^ ndb)
        #xK, yK, zK = approxpolycentroid3(pv, Tri)
        #@show (xK, yK, zK)
        xK, yK, zK = centroids[iel, :]
        #@show (xK, yK, zK)
        hK = diameters[iel]
        x = @view V[:, 1]
        y = @view V[:, 2]
        z = @view V[:, 3] 
        # ------- scaled monomials ----------
        m1(x, y, z) = ones(T, length(x))
        m2(x, y, z) = (x .- xK) / hK
        m3(x, y, z) = (y .- yK) / hK
        m4(x, y, z) = (z .- zK) / hK
        m(x, y, z) = [m1(x, y, z)[1], m2(x, y, z), m3(x, y, z), m4(x , y, z)]
        mc = [m1, m2, m3, m4]
        gradmMat = [0 0 0; 1/hK 0 0; 0 1/hK 0; 0 0 1/hK]
        #-------- transition matrix ----------
        D = hcat(m1(x, y, z), m2(x, y, z), m3(x, y, z), m4(x , y, z))
        # ----------- elliptic projection -------------
        B = zeros(4, Ndof)
        #for s = 1:length(elemf)
        for s = 1:length(elemf)
            # --- information of current face
            # vertices of face
            facesloc = elemf[s]
            P = @view pv[facesloc, :]
            # elliptic projection on the face
            idFace = indexFace[s]   
            Pifs = faceProj[idFace] # the order may be not correct
            # strange sorting.. (not understood completely)
            funik = sort(unique(facesloc))
            idx = [findfirst(isequal(ff), funik) for ff ∈ facesloc]
            Pifs = Pifs[:, idx]
            # normal vector           
            e1 = P[2, :] - P[1, :]
            en = P[1, :] - P[end, :]
            nf = cross(e1, en)
            nf ./= norm(nf)
            # area
            areaf = polyarea(P)
            # --- integral of Pifs
            intFace = [areaf 0 0] * Pifs  # local
            #display(Pifs)
            intProj = zeros(1 , Ndof)     # global adjustment
            faceLocal = elemfLocal[s]      
            intProj[1, faceLocal] .= intFace[:]
            # add grad(m) * nf
            Bf = sum(gradmMat .* nf', dims = 2) * intProj
            #display(intProj)
            B = B + Bf
        end
        # constraint
        Bs = copy(B)  
        Bs[1, :] .= 1 / Ndof
        # consistency relation
        G = B * D
        Gs = Bs * D
        # --------- L2 projection -----------   
        H = zeros(4, 4)
        for i = 1:4
            fun(X) = mc[i](X[1], X[2], X[3])[1] * m(X[1], X[2], X[3])
            H[i, :] = integralPolyhedron(fun, 3, pv, elemf, centroids[iel, :]) 
            ## centroids??
        end
        # --------- local stiffness matrix ---------
        Pis = Gs \ Bs  
        Pi  = D * Pis
        I = Matrix(LinearAlgebra.I, size(Pi, 1), size(Pi, 2))
        #cc = 1.
        #AK  = Pis' * G * Pis + cc * Pis' * H * Pis +
        #      hK * (1 + cc * hK ^ 2) * (I - Pi)' * (I - Pi) 
        AK  = Pis' * G * Pis + hK * (I - Pi)' * (I - Pi) 
        #C = H * (Gs \ Bs)
        #AM  = C' * (H \ C) + hK ^ 3 * (I - Pi)' * (I - Pi) 
        AM  = Pis' * H * Pis + hK ^ 3 * (I - Pi)' * (I - Pi) 
        # --------- assembly index for ellptic projection -----------
        ii[(ia + 1):(ia + Ndof ^ 2)] = repeat(indexDof', Ndof, 1)
        jj[(ia + 1):(ia + Ndof ^ 2)] = repeat(indexDof, Ndof)
        SK[(ia + 1):(ia + Ndof ^ 2)] = AK
        SM[(ia + 1):(ia + Ndof ^ 2)] = AM
        ia += Ndof ^ 2
        # --------- matrix for L2 and H1 error evaluation  ---------
        push!(Ph, copy(Pis))
        push!(elem2dof, copy(indexDof))
    end
    #=
    # octave comparison
     S = SK + SM
    ijs = readdlm("/tmp/kk.mat")[6:end, 1:3]
    @show norm(ijs[:, 1] - ii)
    @show norm(ijs[:, 2] - jj)
    @show norm(ijs[:, 3] - S)
    =#
    IK, JK, IM, JM = ii, jj, ii, jj
    return IK, JK, SK, IM, JM, SM
end
#-------------------------------------------------------------------------------
"""
    vem_eigs3D(filename::String =  "Lpolmesh_coarse", 
               nc::Int64 = 1_00; resolution::Int64 = 400)
               → u, p, t, meshboundary 

N.B. 
* nc = 100, 1_000 or 10_000
* adapted code from the iVEM-1.0 Matlab Toolbox
"""
function vem_eigs3D(filename::String = "Lpolmesh", nc::Int64 = 1_00;
                    resolution::Int64 = 400, visible::Bool = false,
                    numeig::Int64 = 2)
    # computes the virtual element solution of the Poisson problem on 
    # the specified meshprintln(norm.(centers .- c.centers))

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
function solve_eigs3D(IK, JK, SK, IM, JM, SM, pv, meshboundary,
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
function dirichlet3D(points)
    return 0.
end
