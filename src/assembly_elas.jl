function boundary_condition_L_elas(c)
    g = [0.0, 0.0]
    return g
end
function rhs_L_elas(c)
    f = [0.0, 0.0]
    if c[1] > 0.9
        f[1] = -1.0
    end
    return f
end
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-------------------------
function characElt(verts) 
    n, n_sides = size(verts, 1), size(verts, 1)
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
    return area, centroid, diameter
end
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-------------------------
"""
   IK, JK, SK, IF, SF = assembKM_vemKsource_elas(pv, cellsb, rhs, 
                                                 boundary_condition,
                                                 meshboundary)

Assemble stiffness and source terms of virtual element approach
"""
function assembKM_vemKsource_elas(pv, cellsb, rhs, boundary_condition,
                                  meshboundary)
    # material parameters
    lm = 0.5769
    mu = 0.3846
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
    # a utility function for wrapping around a vector
    mod_wrap(x, a) = mod(x - 1, a) + 1 
    # (3*3) matrix D of elasticity over constant strain fields
    Di = [ (2.0 * mu + lm)              lm         0.0;
                       lm   (2.0 * mu + lm)        0.0;
                      0.0              0.0    4.0 * mu]
    for el_id = 1:length(cellsb)
        #println(cellsb[el_id])
        v = cellsb[el_id]      # global indices of the vertices of el 
        vc = pv[v, :]   # na * 2 array of vertices coordinates
        na = length(v)  
        # start computing the geometric information
        vol, bar, diameter = characElt(vc) 
        D = vol * Di
        # basis functions for R = [  (1,0)  ,  (0,1)  , (-(x_1-bar1,x_0-bar0)) ]
        # basis functions for C = [ (x_0-bar0,0) , (0,x_1-bar1) ,(x_1-bar1, 
        #                            x_0-bar0)  ]
        # (2n*3) Matrix WC of components of basis functions over constant  
        # strain fields
        WC = zeros(typein, 2 * na, 3)
        # (2n*3) Matrix WR of components of basis functions over rigid-body 
        # motions
        WR = zeros(typein, 2 * na, 3)
        # (2n*3) matrix NC: (NC)_{2i,l} = (c_k(x_i))_0, (NC)_{2i+1,l} = 
        # (c_k(x_i))_1
        NC = zeros(typein, 2 * na, 3)
        # (2n*3) matrix NR: (NR)_{2i,l} = (r_k(x_i))_0, (NR)_{2i+1,l} = 
        # (r_k(x_i))_1
        NR = zeros(typein, 2 * na, 3)
        for i = 1:na
            # coordinates of vertices i-1, i, i+1
            ci = vc[i, :] # this vertex and its neighbours
            cim = vc[mod_wrap(i - 1, na), :]
            cip = vc[mod_wrap(i + 1, na), :]
            # normal vector at vertex (multiplied by its length)
            nor = [cip[2] - cim[2], - (cip[1] - cim[1])]
         
            q1 = 0.5 * nor[1] / vol
            q2 = 0.5 * nor[2] / vol
            
            WC[2 * i - 1, 1]   = q1
            WC[2 * i - 1, 2]   = 0.0
            WC[2 * i - 1, 3]   = 0.5 * q2
            
            WC[2 * i, 1] = 0.0
            WC[2 * i, 2] = q2
            WC[2 * i, 3] = 0.5 * q1
            
            WR[2 * i - 1, 1]   = 1.0 / na
            WR[2 * i - 1, 2]   = 0.0
            WR[2 * i - 1, 3]   = 0.5 * q2
            
            WR[2 * i, 1] = 0.0
            WR[2 * i, 2] = 1.0 / na
            WR[2 * i, 3] = -0.5 * q1
            
            NC[2 * i - 1, 1] = ci[1] - bar[1]
            NC[2 * i - 1, 2] = 0.0
            NC[2 * i - 1, 3] = ci[2] - bar[2]
            
            NC[2 * i, 1] = 0.0
            NC[2 * i, 2] = ci[2] - bar[2]
            NC[2 * i, 3] = ci[1] - bar[1]
            
            NR[2 * i - 1, 1] = 1.0
            NR[2 * i - 1, 2] = 0.0
            NR[2 * i - 1, 3] = ci[2] - bar[2]
            
            NR[2 * i, 1] = 0.0
            NR[2 * i, 2] = 1.0
            NR[2 * i, 3] = - (ci[1] - bar[1])
        end
        # projection matrix over polynomial displacements
        PP = WC * NC' + WR * NR'
        # projection matrix over higher-order functions
        QP = eye(2 * na) - PP
        # polynomial block in the local stiffness matrix
        PROJ = WC * D * WC'
        # higher-order functions block in the local stiffness matrix
        alphaE = tr(D) / tr(NC' * NC)
        STAB = alphaE * QP * QP'
        # local stiffness matrix
        Ke = PROJ + STAB
        # copy local to global
        # K[vert_ids, vert_ids] += local_stiffness 
        # F[vert_ids] .+= rhs(centroid) * area / n_sides
        for il = 1:na
            ig = v[il]
            for jl = 1:na
                jg = v[jl]
                push!(IK, 2 * ig - 1)
                push!(IK, 2 * ig - 1)
                push!(IK, 2 * ig)
                push!(IK, 2 * ig)

                push!(JK, 2 * jg - 1)
                push!(JK, 2 * jg)
                push!(JK, 2 * jg - 1)
                push!(JK, 2 * jg)

                push!(SK, Ke[2 * il - 1, 2 * jl - 1])
                push!(SK, Ke[2 * il - 1, 2 * jl])
                push!(SK, Ke[2 * il, 2 * jl - 1])
                push!(SK, Ke[2 * il, 2 * jl])
            end
        end
        # update force vector
        f = rhs(bar)
        for il = 1:na
            ig = v[il]
            push!(IF, 2 * ig - 1)
            push!(IF, 2 * ig)
            push!(SF, vol / na * f[1])
            push!(SF, vol / na * f[2])            
        end
    end
    # add boundary contribution in the source term 
    nb = length(meshboundary)
    boundary_vals = [boundary_condition(pv[meshboundary[k], :]) for k = 1:nb]
    # boundary dofs
    Vmeshboundary = vcat(2 * meshboundary .- 1, 2 * meshboundary)
    Jinmesh = indexin(JK, Vmeshboundary)
    Jboundary = findall(Jinmesh .!= nothing)
    for l ∈ Jboundary
        push!(IF, IK[l])
        if mod(Jinmesh[l], 2) == 1
            push!(SF, - SK[l] * boundary_vals[(Jinmesh[l] + 1) ÷ 2][1])
        else
            push!(SF, - SK[l] * boundary_vals[Jinmesh[l] ÷ 2][2])
        end
    end
    return IK, JK, SK, IF, SF
end
#-------------------------------------------------------------------------------
"""
    u, p, t, meshboundary = vem_elas(filename::String = "squarepolmesh_coarse",
                                     nc::Int64 = 1_00;
                                     resolution::Int64 = 400)

N.B. 
* nc = 100, 1_000 or 10_000
* adapted code from the article / matlab's code  "The virtual element method in 50 lines of matlab". See https://arxiv.org/pdf/1604.06021.pdf
"""
function vem_elas(filename::String = "Lpolmesh", nc::Int64 = 1_00;
                  resolution::Int64 = 400)
    # computes the virtual element solution of the Poisson problem on 
    # the specified mesh
    # load the mesh  + pv       : vertices of the cells 
    #                + cellsb   : polygonal cells (⚠ orientation is crucial)
    #                + celssbt  : triangulated cells 
    #                + t        : all triangles 
    #                + (pb, tb) : restricted delaunay mesh 
    mesh_filename = "$(@__DIR__)/../test/data/$(filename)_$(nc).jld2"
    println("read file $(mesh_filename)")
    JLD2.@load(mesh_filename, pv, cellsb, cellsbt, t, pb, tb) 
    # boundary points / rhs  
    rhs, boundary_condition, meshboundary = if occursin("square", mesh_filename)
        rhs_sqr_elas, boundary_condition_sqr_elas, []
    else
        EPS = 1e-4
        rhs_L_elas, boundary_condition_L_elas, findall(pv[:, 2] .> (1 - EPS))
    end
    n_dofs, n_polys = size(pv, 1), 3 # method has 1 degree of freedom per vertex
    n_dofs *= 2    
    u = zeros(n_dofs) # degrees of freedom of the virtual element solution
    # call assemble function
    IK, JK, SK, IF, SF = assembKM_vemKsource_elas(pv, cellsb, rhs,
                                                  boundary_condition, 
                                                  meshboundary)
    K = sparse(IK, JK, SK, n_dofs, n_dofs) 
    F = Vector(sparsevec(IF, SF, n_dofs))
    nb = length(meshboundary)
    boundary_vals = [boundary_condition(pv[meshboundary[k], :]) for k = 1:nb]
    Vboundary_vals = [boundary_vals[k][l] for l = 1:2, k = 1:nb][:]
    # vertices which aren’t on the boundary
    Vmeshboundary = vcat(2 * meshboundary .- 1, 2 * meshboundary)
    internal_dofs =  setdiff(1:n_dofs, Vmeshboundary) 
    # solve
    u[internal_dofs] = K[internal_dofs, internal_dofs] \ F[internal_dofs] 
    u[Vmeshboundary] .= Vboundary_vals # set the boundary values
    # plot
    if resolution > 0
        plotsolution(u, pv, cellsb, mesh_filename = mesh_filename,
                     resolution = resolution)
    end
    return u, pv, t, meshboundary
end