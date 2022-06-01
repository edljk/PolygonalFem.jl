using Test, JLD2, LinearAlgebra
#
# gradient computation by parial automatic differentiation and adjoint methods 
# of  
#           λ(Ω) 
# where 
# K(v) U(v) = λ(Ω) M(v) U(v) the discrete version of - Δu = λ(Ω) u (Dirichlet) 
# and v are the vertices of the voronoi cells 
# 

# load polygonal mesh 
filename, nc = "Lpolmesh", 1_000
mesh_filename = "$(@__DIR__)/data/$(filename)_$(nc).jld2"
println("read file $(mesh_filename)")
JLD2.@load(mesh_filename, pv, bbelem) 
cellsb = bbelem
# boundary points 
fD(x) = 0.
numeig = 1
boundary_condition, meshboundary = if occursin("square", mesh_filename)
    EPS = 1e-4
    fD, union(findall(pv[:, 2] .> (1 - EPS)),
              findall(pv[:, 2] .< EPS),
              findall(pv[:, 1] .> (1 - EPS)),
              findall(pv[:, 1] .< EPS))
else
    EPS = 1e-4
    fD, union(findall(pv[:, 2] .> (1 - EPS)),
              findall(pv[:, 2] .< - 1 + EPS),
              intersect(findall(abs.(pv[:, 2]) .< EPS),
                        findall(pv[:, 1] .>= 0.)),
              findall(pv[:, 1] .> (1 - EPS)),
              findall(pv[:, 1] .< - 1 + EPS),
              intersect(findall(abs.(pv[:, 1]) .< EPS),
                        findall(pv[:, 2] .>= 0.)))
end
np, ε = size(pv, 1), 1e-4
x, dp = pv, 2 * rand(size(pv)...) .- 1
# IMPORTANT Do not move boundary nodes
#dp[meshboundary, :] .= 0.
xp, xm = x .+ ε * dp, x .- ε * dp
# finite differences
f = x -> PolygonalFem.costeigs2D(x, cellsb, boundary_condition, meshboundary,
                                 numeig)[1]
fp, fm = f(xp), f(xm)
@time g = PolygonalFem.∇costeigs2D(x, cellsb, boundary_condition, meshboundary, 
                                   numeig)
@time g = PolygonalFem.∇costeigs2D(x, cellsb, boundary_condition, meshboundary,
                                   numeig)
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