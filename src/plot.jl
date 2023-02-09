function plotunicode_points(p)
    show(UnicodePlots.scatterplot(p[:, 1], p[:, 2], 
                                  marker = :diamond))
    return
end
function plotunicode(U;
                     colormap::Symbol = :viridis, title::String = "")
    np = length(U)
    snp = ceil(Int64, sqrt(np))
    u = reshape(U, snp, snp)
    show(UnicodePlots.heatmap(real.(u), xfact = 2 / (snp - 1), 
                              yfact = 2 / (snp - 1),
                              xoffset = -1., yoffset = -1., 
                              colormap = colormap, title = title,
                              #, colormap = :inferno,
                              width = 30, height = 30))
    print("\n")
    return nothing 
end
#-------------------------------------------------------------------------------
"""
   
   plotmesh(filename::String = "Lpolmesh", nc::Int64 = 1_00;
            dimplot::Bool = true, alphac::Float64 = 0.9,  
            visible::Bool = false, cut::Int64 = 0)
"""
function plotmesh(filename::String = "Lpolmesh", nc::Int64 = 1_00;
                  dimplot::Bool = true, alphac::Float64 = 0.9,  
                  visible::Bool = true, cut::Int64 = 0,
                  randomcolors::Bool = !dimplot)
    mesh_filename = "$(@__DIR__)/../test/data/$(filename)_$(nc).jld2"
    println("read file $(mesh_filename)")
    JLD2.@load(mesh_filename, ps, pv, elem, bbelem, Iv, pb, tb) 
    u = if dimplot 
        pv'[:]
    else 
        pv[:, end]
    end
    dim = size(ps, 2)
    cells = dim == 2 ? bbelem : elem
    bcells = bbelem
    if !visible 
        Iv = 1:size(ps, 1)
    end
    if cut > 0
        Iv = intersect(Iv, findall((ps .- mean(pv, dims = 1))[:, cut] .> 0))
    end
    plotsolution(u, pv, cells, bcells, Iv, ps, wireframe = true, 
                 visible = visible,
                 alphac = alphac, randomcolors = randomcolors)
    nothing
end
#-------------------------------------------------------------------------------
function plotsolution3D(uin, pv, cells, bcells, Iv, ps;
                        wireframe::Bool = false, randomcolors::Bool = false,
                        alphac::Float64 = 0.9, visible::Bool = false,
                        resolution::Int64 = 400, mesh_filename::String = "")
    np, m = length(cells), GeometryBasics.Mesh[] 
    GLMakie.destroy!(GLMakie.Screen())
    colS = range(Colors.HSV(0, 1, 1), stop = Colors.HSV(330, 1, 1), 
                 length = 64)
    colmapS = [convert(Colors.RGB{Float32}, colS[k]) for  k = 1:length(colS)] 
    dim = length(uin) ÷ size(pv, 1)
    fig = GLMakie.Figure(resolution = (dim * resolution, resolution)) 
    println("build meshes in plotsolution3D")
    pts_wire = Vector{Vector{Float64}}()
    for k ∈ Iv #1:np
        Ik = unique(cells[k][:])
        pk = pv[Ik, :]
        t = convexhull(pk)[1]
        if wireframe
            fG = bcells[k]
            #fG = convexhull_3Dfaces(pk)[1]            
            for l = 1:length(fG)
                for m = 1:length(fG[l])
                    if m == 1
                        push!(pts_wire, pk[fG[l][1], :])
                        push!(pts_wire, pk[fG[l][end], :])
                    else
                        push!(pts_wire, pk[fG[l][m - 1], :])
                        push!(pts_wire, pk[fG[l][m], :])
                    end
                end
            end
        end
        tG = [TriangleFace{Int32}(t[l, 1], t[l, 2], t[l, 3]) 
                                                          for l = 1:size(t, 1)]
        pG =  [Point3f0(pk[l, 1], pk[l, 2], pk[l, 3]) for l = 1:size(pk, 1)]
        push!(m, normal_mesh(copy(pG), copy(tG)))
    end
    if wireframe  # add edges of polygons 
        npw = length(pts_wire)
        x = [pts_wire[k][1] for k = 1:npw]
        y = [pts_wire[k][2] for k = 1:npw]
        z = [pts_wire[k][3] for k = 1:npw]
    end
    # plot 
    axlist = Any[]
    for k = 1:dim 
        ax = L = GLMakie.Axis3(fig[1, 2 * k - 1])
        push!(axlist, ax)
        u, title = if dim == 1
            uin, ""
        elseif k == 1
            uin[1:3:end], "ux"
        elseif k == 2
            uin[2:3:end], "uy"
        elseif k == 3
            uin[3:3:end], "uz"
        end
        if length(title) > 0 
            ax.title = title
        end
        xu = range(minimum(u), stop = maximum(u), length = length(colS))
        il = LinearInterpolation(xu, colmapS)
        colors = randomcolors ? [RGBA(rand(3)..., alphac) for _ ∈ Iv] : [il(mean(u[cells[l]])) for l ∈ Iv]
        Makie.mesh!(L, m, color = colors, shading = false)
        #Makie.scatter!(L, p[Ipv, 1], p[Ipv, 2], p[Ipv, 3], markersize = 5_000)
        #Makie.scatter!(L, p[Iv, 1], p[Iv, 2], p[Iv, 3], markersize = 5_000)
        ax.aspect = :data
        if wireframe  # add edges of polygons 
            Makie.linesegments!(L, x, y, z, color = :black)
        end
        Colorbar(fig[1, 2 * k], limits = (minimum(u), maximum(u)), 
                 colormap = :hsv,
                 size = 10,  height = Relative(3/4))
    end
    if dim > 1
    #    Main.Makie.linkaxes!(axlist...)
    end
    display(fig)
    return
end
#-------------------------------------------------------------------------------
function plotsolution(uin, pv, cells, bcells, Iv, ps;
                      alphac::Float64 = 0.9, visible::Bool = false,
                      resolution::Int64 = 400,  wireframe::Bool = false,
                      randomcolors::Bool = false,
                      mesh_filename::String = "")
    if size(pv, 2) == 3
        if !scipy
            error("3D plots only available if scipy is installed in conda..")
        end
        plotsolution3D(uin, pv, cells, bcells, Iv, ps, resolution = resolution, 
                       alphac = alphac, mesh_filename = mesh_filename, 
                       wireframe = wireframe, visible = visible,
                       randomcolors = randomcolors)
        return
    end
    GLMakie.destroy!(GLMakie.Screen())
    dim = length(uin) ÷ size(pv, 1)
    fig = GLMakie.Figure(resolution = (dim * resolution, resolution)) 
    P = Polygon[]
    for k ∈ Iv
        cc = [cells[k]..., cells[k][1]]
        push!(P, Polygon([Point2(pv[m, 1], 
                                 pv[m, 2]) for m ∈ cc])) # u[m]
    end
    colS = range(Colors.HSV(0, 1, 1), stop = Colors.HSV(330, 1, 1), 
                 length = 64)
    colmapS = [convert(Colors.RGB{Float32}, colS[k]) for  k = 1:length(colS)]
    dim = length(uin) ÷ size(pv, 1)
    axlist = Any[]
    for k = 1:dim 
        ax = L = GLMakie.Axis(fig[1, 2 * k - 1])
        push!(axlist, ax)
        ax.aspect = GLMakie.DataAspect()
        u, title = if dim == 1
            uin, ""
        elseif k == 1
            uin[1:2:end], "ux"
        elseif k == 2
            uin[2:2:end], "uy"
        end
        if length(title) > 0 
            ax.title = title
        end
        xu = range(minimum(u), stop = maximum(u), length = length(colS))
        il = LinearInterpolation(xu, colmapS)
        colors = randomcolors ? [RGBA(rand(3)..., alphac) for _ ∈ Iv] : [il(mean(u[cells[l]])) for l ∈ Iv]
        Makie.poly!(ax, P, 
                    strokecolor = :black, strokewidth = 2., 
                    color = colors) 
        Colorbar(fig[1, 2 * k], limits = (minimum(u), maximum(u)), 
                 colormap = :hsv,
                 size = 10,  height = Relative(3/4))
        #Makie.scatter!(L, ps[Iv, 1], ps[Iv, 2], markersize = 10)
        if occursin("square", mesh_filename)
            xlims!(ax, (0, 1.)) 
            ylims!(ax, (0, 1.))
        elseif occursin("L", mesh_filename)
            xlims!(ax, (-1, 1.)) 
            ylims!(ax, (-1, 1.))
        end
    end
    if dim > 1
        Makie.linkaxes!(axlist...)
    end
    display(Makie.current_figure())
    return
end
