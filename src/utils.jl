const pyspatial = PyNULL()
function __init__()
    copy!(pyspatial,  pyimport("scipy.spatial"))
end
# is scipy available?
const scipy = pyspatial == PyNULL() ? false : true
#-------------------------------------------------------------------------------
function edgestoloop(el::Array{Int64, 2})
    elcol = reshape(el, length(el), 1)
    listep, I, J = uniquerows(elcol)
    J = reshape(J, size(el))
    ptoe = zeros(Int64, length(I), 2)
    for k = 1:size(el, 1)
        for l = 1:2
            if ptoe[J[k, l], 1] == 0
                ptoe[J[k, l], 1] = k
            else
                ptoe[J[k, l], 2] = k
            end
        end
    end
    if minimum(ptoe) == 0
        error("Pb I in edgestoloop")
    end
    loop = zeros(Int64, length(I)); loop[1:2] = J[1, :]; ec = 1 # courant edge
    for k = 3:length(I)
        ec = ec == ptoe[loop[k - 1], 1] ? ptoe[loop[k - 1], 2] : ptoe[loop[k - 1], 1]
        #ec = setdiff(ptoe[loop[k - 1], :], ec) # actualize courant edge
        # setdiff(J[ec, :], loop[k - 1])[1]
        loop[k] = J[ec, 1] == loop[k - 1] ? J[ec, 2] : J[ec, 1]
    end
    if minimum(loop) == 0
        error("Pb II in edgestoloop")
    end
    return elcol[I[loop]]
end
#-------------------------------------------------------------------------------
"""
    Tri, pyhull = convexhull(p::Array{Float64, 2}, 
                             qhull_option::String = "QJ Pp")
    convexhull()
Call spatial scipy library to compute convehull with qhull.

WARNING : there is perhaps a problem of orientation of triangles!
"""
function convexhull(p::Array{Float64, 2}, qhull_option::String = "QJ Pp")
    pyhull = pyspatial.ConvexHull(p, qhull_options = qhull_option)
    Tri = convert(Array{Int64, 2}, pyhull.simplices) .+ 1
    return Tri, pyhull
end
#-------------------------------------------------------------------------------
function convexhull_3Dfaces(p::Array{Float64, 2}, 
                            preci::Int64 = 6,
                            qhull_option::String = "QJ Pp")
    if size(p, 2) != 3
        error("convexhull_3Dfaces relevant only in dimension 3..")
    end
    Tri, pyhull = convexhull(p, qhull_option) 
    Tri = convert(Array{Int64, 2}, pyhull.simplices) .+ 1
    eq = pyhull.equations # triangular faces equations
    _, Iu, Ju = uniquerows(eq, preci)
    # reconstruct faces 
    faces = Vector{Vector{Int64}}()
    M = zeros(3, 3)
    for k ∈ unique(Ju)
        tk = findall(Ju .== k)
        if length(tk) == 1 # triangular face
            push!(faces, Tri[tk[1], :])
        else # polygonal face
            n = eq[tk[1], 1:3]
            # project to planar configuration
            M[:, 3] .= n 
            if maximum(abs.(n[1:2])) > 1e-8
                M[:, 2] .= [- n[2], n[1], 0.]
            else
                M[:, 2] .= [0., - n[3], n[2]]
            end
            M[:, 1] .= cross(M[:, 2], M[:, 3])
            Ik = unique(Tri[tk, :][:])
            Minv = inv(M)
            p2D = copy((Minv * p[Ik, :]')'[:, 1:2])
            push!(faces, Ik[Main.edgestoloop(convexhull(p2D)[1])])
        end
    end
    return faces, pyhull
end
#-------------------------------------------------------------------------------
"""
    Au, Iu, Ju = uniquerows(A::AbstractArray{T}, preci::Int64 = 6)
Matlab like unique_rows function

N.B. 
* use preci = 0 for integers
* maximum(abs.((A[Iu, :])[Ju, :] .- A)) == 0
"""
function uniquerows(A::AbstractArray{T}, preci::Int64 = 0) where T
    Ar = if preci > 0 # use preci = 0 for integers
        round.(Int64, A * 10 ^ preci) 
    else
        A
    end
    Ju = groupslices(Ar, dims = 1)
    p = sortperm(Ju)
    Iu, c = [Ju[p[1]],], 1
    Jun = similar(Ju)
    Jun[p[1]] = c
    for k = 2:length(Ju)
        if Ju[p[k]] > Ju[p[k - 1]]
            push!(Iu, Ju[p[k]])
            c += 1
        end
        Jun[p[k]] = c 
    end
    #@test maximum(abs.((Ar[Iu, :])[Jun, :] .- Ar)) == 0
    return A[Iu, :], Iu, Jun
end