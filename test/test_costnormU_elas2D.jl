using Test, JLD2, LinearAlgebra
#
# gradient computation by parial automatic differentiation and adjoint methods 
# of 
#           E(v) = ∫ |U(v)|² dx 
# where 
# A(v) U(v) = F(v) the discrete version of - div(M ∇ Ua) = Fa
# and v are the vertices of the voronoi cells 
# 

# load polygonal mesh 
filename, nc = "Lpolmesh_coarse", 100
mesh_filename = "$(@__DIR__)/data/$(filename)_$(nc).jld2"
println("read file $(mesh_filename)")
JLD2.@load(mesh_filename, pv, cellsb, cellsbt, t, pb, tb) 
# boundary points / rhs 
rhs, boundary_condition, meshboundary = if occursin("square", mesh_filename)
    PolygonalFem.rhs_sqr, PolygonalFem.boundary_condition_sqr, []
else
    EPS = 1e-4
    PolygonalFem.rhs_L_elas, PolygonalFem.boundary_condition_L_elas, 
    findall(pv[:, 2] .> (1 - EPS))
end
np, ε = size(pv, 1), 1e-4
x, dp = pv, 2 * rand(size(pv)...) .- 1
# IMPORTANT Do not move boundary nodes
dp[meshboundary, :] .= 0.
xp, xm = x .+ ε * dp, x .- ε * dp
# finite differences
f = x -> PolygonalFem.costnormU_elas2D(x, cellsb, t, rhs, boundary_condition, 
                                       meshboundary)
fp, fm = f(xp), f(xm)
@time g = PolygonalFem.∇costnormU_elas2D(x, cellsb, t, rhs, boundary_condition, 
                                         meshboundary)
@time g = PolygonalFem.∇costnormU_elas2D(x, cellsb, t, rhs, boundary_condition,
                                         meshboundary)
println("")
println(fp)
println(fm)
println(" ")
println((fp - fm) / (2 * ε))
println(dot(g, dp))
println("")
@testset "test eval gradient dot function for polygonal costnormU_elas2D" begin 
    @test abs(dot(g, dp) - (fp - fm) / (2 * ε)) < 1e-4
end