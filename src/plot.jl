
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
function plotmesh(filename::String = "Lpolmesh", nc::Int64 = 1_00;
                  dimplot::Bool = true)
    mesh_filename = "$(@__DIR__)/../test/data/$(filename)_$(nc).jld2"
    println("read file $(mesh_filename)")
    JLD2.@load(mesh_filename, pv, cellsb, cellsbt, t, pb, tb) 
    u = if dimplot
        pv'[:]
    else 
        pv[:, end]
    end
    plotsolution(u, pv, cellsb, wireframe = true)
    nothing
end
#-------------------------------------------------------------------------------
function plotsolution3D(uin, p, cellsb;
                        wireframe::Bool = false,
                        alphac::Float64 = 0.9, 
                        resolution::Int64 = 400, mesh_filename::String = "")
    np, m = length(cellsb), GeometryBasics.Mesh[] 
    GLMakie.destroy!(GLMakie.global_gl_screen())
    colS = range(Colors.HSV(0, 1, 1), stop = Colors.HSV(330, 1, 1), 
                 length = 64)
    colmapS = [convert(Colors.RGB{Float32}, colS[k]) for  k = 1:length(colS)] 
    dim = length(uin) ÷ size(p, 1)
    fig = GLMakie.Figure(resolution = (dim * resolution, resolution)) 
    println("build meshes in plotsolution3D")
    pts_wire = Vector{Vector{Float64}}()
    for k = 1:np
        Ik = unique(cellsb[k][:])
        pk = p[Ik, :]
        t = Main.convexhull(pk)[1]
        if scipy 
            fG = convexhull_3Dfaces(pk)[1]
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
    if wireframe && scipy # add edges of polygons 
        npw = length(pts_wire)
        x = [pts_wire[k][1] for k = 1:npw]
        y = [pts_wire[k][2] for k = 1:npw]
        z = [pts_wire[k][3] for k = 1:npw]
    end

    # plot 
    for k = 1:dim 
        ax = L = GLMakie.Axis3(fig[1, 2 * k - 1])
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
        colors = false ? [RGBA(rand(3)..., alphac) for _ = 1:np] : [il(mean(u[cellsb[l]])) for l = 1:np]
        Makie.mesh!(L, m, color = colors, shading = false)
        ax.aspect = :data
        if wireframe && scipy # add edges of polygons 
            Makie.linesegments!(L, x, y, z, color = :black)
        end
        Colorbar(fig[1, 2 * k], limits = (minimum(u), maximum(u)), 
                 colormap = :hsv,
                 size = 10,  height = Relative(3/4))
    end
    display(fig)
    return
end
#-------------------------------------------------------------------------------
function plotsolution(uin, p, cellsb;
                      resolution::Int64 = 400,  wireframe::Bool = false,
                      mesh_filename::String = "")
    if size(p, 2) == 3
        plotsolution3D(uin, p, cellsb, resolution = resolution, 
                       mesh_filename = mesh_filename, wireframe = wireframe)
        return
    end
    GLMakie.destroy!(GLMakie.global_gl_screen())
    dim = length(uin) ÷ size(p, 1)
    fig = GLMakie.Figure(resolution = (dim * resolution, resolution)) 
    P = Polygon[]
    for k = 1:length(cellsb)
        cc = [cellsb[k]..., cellsb[k][1]]
        push!(P, Polygon([Point2(p[m, 1], 
                                 p[m, 2]) for m ∈ cc])) # u[m]
    end
    colS = range(Colors.HSV(0, 1, 1), stop = Colors.HSV(330, 1, 1), 
                 length = 64)
    colmapS = [convert(Colors.RGB{Float32}, colS[k]) for  k = 1:length(colS)]
    dim = length(uin) ÷ size(p, 1)
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
        Makie.poly!(ax, P, 
                    strokecolor = :black, strokewidth = 2., 
                    color = [il(mean(u[cellsb[k]])) for k = 1:length(P)])
    
        Colorbar(fig[1, 2 * k], limits = (minimum(u), maximum(u)), 
                 colormap = :hsv,
                 size = 10,  height = Relative(3/4))
        if occursin("square", mesh_filename)
            xlims!(ax, (0, 1.)) 
            ylims!(ax, (0, 1.))
        elseif occursin("L", mesh_filename)
            xlims!(ax, (-1, 1.)) 
            ylims!(ax, (-1, 1.))
        end
    end
    if dim > 1
        Main.Makie.linkaxes!(axlist...)
    end
    display(Makie.current_figure())
    return
end
