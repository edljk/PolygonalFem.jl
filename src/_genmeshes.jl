"""
   FOR PRE PROCESSING ONLY (generates files in test/data)
"""
function _genpolmeshes(; numb::Int64 = 2, 
                         nbit::Union{Vector{Int64}, Int64} = [100,],
                         np::Union{Vector{Int64}, Int64} = [100, 1_000, 10_000],
                         drawvoronoi::Bool = false)
    support = if numb == 2
        "square"
    elseif numb == 4
        "flower"
    elseif numb == 9
        "L"
    end  
    # gen regular polygonal meshes using Lloyd
    for itmax ∈ nbit
        strit = if itmax == 1 
            "_coarse" 
        elseif itmax < 1_000 
            "_lesscoarse"
        else
             ""
        end
        for npl ∈ np
            p, pv, cellsb, cellsbt, t, pc, pb, tb = Main.ConvexTools.lloyd_geogram(npl, numb = numb, nbit = itmax, 
                          drawvoronoi = drawvoronoi)
            # save 
            mesh_filename = "$(@__DIR__)/../test/data/$(support)polmesh$(strit)_$(npl).jld2"
            println(mesh_filename)
            JLD2.save(mesh_filename, "pv", pv, "cellsb", cellsb, 
                      "cellsbt", cellsbt, "pb", pb, "tb", tb, "t", t)
        end
    end
    nothing
end
#-------------------------------------------------------------------------------
"""
   FOR PRE PROCESSING ONLY (generates files in test/data)

CALL: _genmeshes(2, np = [100, 1_000, 10_000])
"""
function _genmeshes(; numb::Int64 = 2, 
                      np::Union{Vector{Int64}, Int64} = [100, 1_000, 10_000],
                      drawvoronoi::Bool = false)
    support = if numb == 2
        "square"
    elseif numb == 4
        "flower"
    elseif numb == 9
        "L"
    end  
    # gen regular meshes calling triangle
    for npl ∈ np
        _, ptri, tri = Main.MeshTools.triangle(numb, nbtri = 2 * npl)
        if numb == 2 # prefere a unit uncentered square
            ptri .+= 1
            ptri ./= 2
        end
        # keep only interior points for voronoi cells
        Ik = setdiff(1:size(ptri, 1), unique(btri(tri)[:]))
        p, pv, cellsb, cellsbt, t, pc, pb, tb = Main.ConvexTools.lloyd_geogram(npl, numb = numb, nbit = 1, p = ptri[Ik, :], drawvoronoi = drawvoronoi)
        # save 
        mesh_filename = "$(@__DIR__)/../test/data/$(support)mesh_$(npl).jld2"
        println(mesh_filename)
        JLD2.save(mesh_filename, "pv", pv, "cellsb", cellsb, 
                  "cellsbt", cellsbt, "pb", pb, "tb", tb, "t", t, 
                  "ptri", ptri, "tri", tri)
    end
    nothing
end
#-------------------------------------------------------------------------------
"""
   FOR POST PROCESSING ONLY
"""
function _genfigures()
    for file ∈ ["Lpolmesh", "Lpolmesh_lesscoarse", "Lpolmesh_coarse",
                "squarepolmesh", "squarepolmesh_lesscoarse", 
                "squarepolmesh_coarse"]
        for np ∈ [100, 1000, 10000]
            mesh_filename = "$(@__DIR__)/../test/data/$(file)_$(np).jld2"
            vem(file, np, resolution = 1_000)
            figfile = replace(replace(mesh_filename, "jld2" => "png"), 
                                                     "data" => "figures")
            println("save figure to $(figfile)")   
            Makie.save(figfile, Makie.current_figure())
        end
    end
    nothing 
end
#-------------------------------------------------------------------------------
"""
   FOR POST PROCESSING ONLY
"""
function _gensolutions()
    sols = Dict()
    for file ∈ ["Lpolmesh", "Lpolmesh_lesscoarse", "Lpolmesh_coarse",
                "squarepolmesh", "squarepolmesh_lesscoarse", 
                "squarepolmesh_coarse"]     
        for np ∈ [100, 1000, 10000]
            mesh_filename = "$(@__DIR__)/../test/data/$(file)_$(np).jld2"
            u, pv, t, Ib = vem(file, np, resolution = 0)
            sols[mesh_filename] = (u, pv, t, Ib)
            println("solution of $(mesh_filename)")   
        end
    end
    @save "/tmp/sols.jld2" sols 
    nothing 
end
#-------------------------------------------------------------------------------
"""
   FOR POST PROCESSING ONLY
"""
function _generrors()
    @load "/tmp/sols.jld2" 
    fsolsquare = "$(@__DIR__)/../test/data/squaresolution_100000.jld2"
    fsolL = "$(@__DIR__)/../test/data/Lsolution_100000.jld2"
    dsolsquare, dsolL =  FileIO.load(fsolsquare), FileIO.load(fsolL)
    errors = Dict()
    for ff ∈ keys(sols)
        println(ff)
        dsol = if occursin("square", ff)
            dsolsquare
        else
            dsolL
        end
        # compute error in ||.||∞ for inside points 
        Inb = setdiff(1:size(sols[ff][2], 1), sols[ff][4])
      
        Itri, Iclose, baryc = Main.barycentric_coordinates(dsol["ptri"], 
                                                           dsol["tri"], 
                                                           sols[ff][2])
        Inb = findall(Itri .> 0)
        println("percentage of interior points ", length(Inb) / size(Itri, 1) * 100)
        ui = sum(dsol["u"][dsol["tri"][Itri[Inb], k]] .* baryc[Inb, k] for k = 1:3)
        err = maximum(abs.(ui .- sols[ff][1][Inb]))
        err = norm(ui .- sols[ff][1][Inb]) / sqrt(length(Inb))
        nbdof = length(sols[ff][1])
        errors[ff] = (err, nbdof)
    end
    @save "/tmp/errors.jld2" errors 
    nothing
end
#-------------------------------------------------------------------------------
"""
   FOR POST PROCESSING ONLY
"""
function _plotorder()
    @load "/tmp/errors.jld2"
    Main.closeall()
    f, L = Main.figure(1, square_axis = false, resolution = (400, 800))
    plotsl, plotsm = Any[], Any[]
    markers = [:circle, :rect, :utriangle, :dtriangle]
    allfiles = ["Lpolmesh", "Lpolmesh_lesscoarse", "Lpolmesh_coarse"]
    allfiles = ["squarepolmesh", "squarepolmesh_lesscoarse","squarepolmesh_coarse"]  
    cm = 1
    colors = [RGB(0.8, 0.2, 0.1), RGB(0.8, 0.8, 0.1), RGB(0.1, 0.8, 0.2)]
    for file ∈ allfiles
        x, y = Float64[], Float64[] 
        for np ∈ [100, 1000, 10000]
            mesh_filename = "$(@__DIR__)/../test/data/$(file)_$(np).jld2"
            err, nbdof = errors[mesh_filename]
            push!(x, log(nbdof))
            push!(y, log(err))
            
        end
        println(x)
        println(y)
        cr = colors[cm] #rand(RGB)
        push!(plotsl, Makie.lines!(x, y, linewidth = 5, color = cr))
        push!(plotsm, Makie.scatter!(x, y, color = cr, marker = markers[cm],
              markersize = 12))
        cm += 1
    end
    ax = Axis(f[1, 1])
   # Makie.axislegend(ax, [[plotsl[k], plotsm[k]] for k = 1:length(allfiles)],
   # allfiles, "Selected Dots", position = :rb,
   # orientation = :horizontal)
   Makie.Legend(f[2, 1], [[plotsl[k], plotsm[k]] for k = 1:length(allfiles)],
                 allfiles, tellwidth = false, tellheight = true,
                 orientation = :horizontal)
    Main.GraphicTools.finalizescene()
    nothing
end