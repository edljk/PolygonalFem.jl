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
function plotsolution(u, p, cellsb;
                      resolution::Int64 = 400, mesh_filename::String = "")
    GLMakie.destroy!(GLMakie.global_gl_screen())
    fig = GLMakie.Figure(resolution = (resolution, resolution))
    ax = L = GLMakie.Axis(fig[1, 1])
    ax.aspect = GLMakie.DataAspect() 
    P = Polygon[]
    for k = 1:length(cellsb)
        cc = [cellsb[k]..., cellsb[k][1]]
        push!(P, Polygon([Point2(p[m, 1], 
                                 p[m, 2]) for m ∈ cc])) # u[m]
    end
    colS = range(Colors.HSV(0, 1, 1), stop = Colors.HSV(330, 1, 1), 
                 length = 64)
    colmapS = [convert(Colors.RGB{Float32}, colS[k]) for  k = 1:length(colS)]
    xu = range(minimum(u), stop = maximum(u), length = length(colS))
    il = LinearInterpolation(xu, colmapS)
    Makie.poly!(ax, P, 
                strokecolor = :black, strokewidth = 2., 
                color = [il(mean(u[cellsb[k]])) for k = 1:length(P)])
    
    Colorbar(fig[1, 2], limits = (minimum(u), maximum(u)), colormap = :hsv,
              size = 10,  height = Relative(3/4))
    if occursin("square", mesh_filename)
        xlims!(ax, (0, 1.)) 
        ylims!(ax, (0, 1.))
    elseif occursin("L", mesh_filename)
        xlims!(ax, (-1, 1.)) 
        ylims!(ax, (-1, 1.))
    end
    display(Makie.current_figure())
    return
end
