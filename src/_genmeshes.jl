"""
   ONLY FOR PRE PROCESSING (generates files in test/data)
"""
function _genmeshes(; numb::Int64 = 2, 
                      nbit::Union{Vector{Int64}, Int64} = [100,],
                      np::Union{Vector{Int64}, Int64} = [100, 1_000, 10_000])
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
            p, pv, cellsb, cellsbt, t, pc, pb, tb = Main.ConvexTools.lloyd_geogram(npl, numb = numb, nbit = itmax, drawvoronoi = false)
            # save 
            mesh_filename = "$(@__DIR__)/../test/data/$(support)polmesh$(strit)_$(npl).jld2"
            println(mesh_filename)
            JLD2.save(mesh_filename, "pv", pv, "cellsb", cellsb, 
                      "cellsbt", cellsbt, "pb", pb, "tb", tb, "t", t)
        end
    end
    nothing
end
