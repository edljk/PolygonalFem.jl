using Test, JLD2
#
# gradient computation by parial automatic differentiation and adjoint merthods 
# of 
#           E(v) = ∫ |U(v)|² dx 
# where 
# A(v) U(v) = F(v) the discrete version of - Δ u = f
# and v are the vertices of the voronoi cells 
# 

# load polygonal mesh 
filename, nc = "squarepolmesh_coarse", 100
mesh_filename = "$(@__DIR__)/data/$(filename)_$(nc).jld2"
println("read file $(mesh_filename)")
JLD2.@load(mesh_filename, pv, cellsb, cellsbt, t, pb, tb) 
 # boundary points / rhs 
 rhs, boundary_condition = if occursin("square", mesh_filename)
    PolygonalFem.rhs_sqr, PolygonalFem.boundary_condition_sqr
else
    PolygonalFem.rhs_L, PolygonalFem.boundary_condition_L
end

np, ε = size(pv, 1), 1e-3
x, dp = pv, 2 * rand(size(pv)...) .- 1
Ib = unique(PolygonalFem.btri(t)[:])
xp, xm = x .+ ε * dp, x .- ε * dp
# finite differences
f = x -> PolygonalFem.costnormU(x, cellsb, t, rhs, boundary_condition, Ib)
fp, fm = f(xp), f(xm)
@time g = PolygonalFem.∇costnormU(x, cellsb, t, rhs, boundary_condition, Ib)
@time g = PolygonalFem.∇costnormU(x, cellsb, t, rhs, boundary_condition, Ib)
println("")
println(fp)
println(fm)
println(" ")
println((fp - fm) / (2 * ε))
println(dot(g, dp))
println("")
@testset "test eval direct gradient dot function for costnormU" begin 
    @test abs(dot(g, dp) - (fp - fm) / (2 * ε)) < 1e-4
end