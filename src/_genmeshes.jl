"""
   ONLY FOR PRE PROCESSING (generates files in test/data)
"""
function _genmeshes(; numb::Int64 = 2, nbit::Int64 = 1_000)
    support = if numb == 2
        "square"
    elseif numb == 4
        "flower"
    end 
    # gen regular polygonal meshes using Lloyd
    for np ∈ [10,] #[100, 1_000, 10_000]
        p, pv, cellsb, cellsbt, t, pc, pb, tb = Main.ConvexTools.lloyd_geogram(np, numb = numb, nbit = nbit, drawvoronoi = false)
        # save 
        mesh_filename = "$(@__DIR__)/../test/data/$(support)polmesh_$(np).jld2"
        println(mesh_filename)
        JLD2.save(mesh_filename, "pv", pv, "cellsb", cellsb, 
                  "cellsbt", cellsbt, "pb", pb, "tb", tb, "t", t)
    end
    nothing
end
