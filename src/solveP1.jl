function P1(filename::String = "squaremesh", nc::Int64 = 1_00;
            resolution::Int64 = 400)
    # computes the virtual element solution of the Poisson problem on 
    # the specified mesh
    # load the mesh  + ptri     : points triangulation
    #                + tri      : triangulation from `triangle``
    #                + pv       : vertices of the cells 
    #                + cellsb   : polygonal cells (⚠ orientation is crucial)
    #                + celssbt  : triangulated cells 
    #                + t        : all triangles 
    #                + (pb, tb) : restricted delaunay mesh 
    mesh_filename = "$(@__DIR__)/../test/data/$(filename)_$(nc).jld2"
    println("read file $(mesh_filename)")
    @time JLD2.@load(mesh_filename, ptri, tri ) #pv, cellsb, cellsbt, t, pb, tb)
    println("done")
    # boundary codition and source
    rhs, boundary_condition = if occursin("square", mesh_filename)
        rhs_sqr, boundary_condition_sqr
    else
        rhs_L, boundary_condition_L
    end
    # assemble matrices 
    K, M = assembKM_P12D(ptri, tri)
     # Dirichlet conditions
    np, Ib = size(ptri, 1), unique(btri(tri)[:])
    pfact = 1e8
    K += sparse(Ib, Ib, fill(pfact, length(Ib)), np, np)
    F = M * rhs(ptri) 
    vb = boundary_condition(ptri[Ib, :])
    F[Ib] .+= pfact * vb
    u = K \ F 
    println(maximum(abs.(u[Ib] .- vb)))
    # plot 
    Main.closeall()
    f, L = Main.figure(1, fig3D = true)
    Main.plot_t(hcat(ptri, u), tri, representation = "meshsurface", 
                scalars = u, colormap = "hsv", colorw = (0., 0, 0), 
                scene = L[1])
    return u, ptri, tri, Ib 
end