#-------------------------------------------------------------------------------
"""
    centers_momentum_geogram(ps, pv::Matrix{T}, elem) where T
                             → measures, centers, momentum, diameters 

"""
function centers_momentum_geogram(ps, pv::Matrix{T}, elem) where T
    np, dim = length(elem), size(pv, 2)
    # cell measures, centers and momentum of order 2
    measures = [zeros(T, 0) for k = 1:np] 
    centers = [zeros(T, dim) for k = 1:np]
    momentum =  [zeros(T, dim) for k = 1:np]
    diameters = zeros(T, np)
    A, AB, AC, AD = zeros(T, dim), zeros(T, dim), zeros(T, dim), zeros(T, dim)
    pw = [zeros(T, dim) for _ = 1:(dim + 1)]
    a, b = (5 - sqrt(5)) / 20, (5 + 3 * sqrt(5)) / 20
    for k = 1:np
        measures[k] = if !isempty(elem[k]) > 0
            direct_meshmeasures(pv, elem[k])                       
        else
            zero(T)
        end
        diameters[k] = maximum(norm(pv[l, :] - pv[m, :]) for l ∈ elem[k], 
                                                             m ∈ elem[k])
        for l = 1:size(elem[k], 1)
            if dim == 2
                A  .= pv[elem[k][l, 1], :]
                AB .= pv[elem[k][l, 2], :] .- A
                AC .= pv[elem[k][l, 3], :] .- A
                pw[1] .= A .+     AB / 6. .+     AC / 6.
                pw[2] .= A .+ 2 * AB / 3. .+     AC / 6.
                pw[3] .= A .+     AB / 6. .+ 2 * AC / 3.     
                for m = 1:(dim + 1) 
                    centers[k] .+= pw[m] / 3 * measures[k][l]
                    momentum[k] .+= (pw[m] .^ 2) / 3. * measures[k][l]
                end
            elseif dim == 3
                A  .= pv[elem[k][l, 1], :]
                AB .= pv[elem[k][l, 2], :] .- A
                AC .= pv[elem[k][l, 3], :] .- A
                AD .= pv[elem[k][l, 4], :] .- A
                pw[1] .= A .+ a * AB .+ a * AC .+ a * AD
                pw[2] .= A .+ a * AB .+ b * AC .+ a * AD
                pw[3] .= A .+ a * AB .+ a * AC .+ b * AD 
                pw[4] .= A .+ b * AB .+ a * AC .+ a * AD 
                for m = 1:(dim + 1)       
                    centers[k]  .+= pw[m] / 4. * measures[k][l]
                    momentum[k] .+= (pw[m] .^ 2) / 4. * measures[k][l]
                end
            end
        end
        centers[k] ./= sum(measures[k])
    end
    for k = 1:np 
        momentum[k] .-= sum(measures[k]) * (2 * ps[k, :] .* centers[k] .- 
                                             ps[k, :] .^ 2)
    end
    return sum.(measures), [centers[k][l] for k = 1:length(centers), l = 1:dim],
           momentum, diameters
end
#-------------------------------------------------------------------------------
"""
    direct_meshmeasures(p::Matrix{T}, elem) where T →  meas

Function to compute edges, triangles, tetrahedrons measures. May be used for ForwardDiff-erentiation
"""
function direct_meshmeasures(p::Matrix{T}, elem) where T
    np, dim = size(p)
    nelem, dime = size(elem, 1), size(elem, 2)
    meas = zeros(T, nelem)
    for k = 1:nelem
        if dime == 2     # edges
            meas[k] = norm(p[elem[k, 1], :] - p[elem[k, 2], :]) 
        elseif dime == 3 # triangles
            A  = p[elem[k, 1], :] 
            AB = p[elem[k, 2], :]  .- A
            AC = p[elem[k, 3], :]  .- A
            meas[k] = if dim == 2
                abs(AB[1] * AC[2] - AB[2] * AC[1]) / 2.
            elseif dim == 3
                norm(cross(AB, AC)) / 2.
            end
        elseif dime == 4  # tetrahedrons
            d12 = p[elem[k, 2], :] - p[elem[k, 1], :]
            d13 = p[elem[k, 3], :] - p[elem[k, 1], :]
            d14 = p[elem[k, 4], :] - p[elem[k, 1], :]
            meas[k] = abs(sum(cross(d12, d13) .* d14)) / 6.
        end
    end
    return meas 
end
#-------------------------------------------------------------------------------
function faceEllipticProjection(P)
    # compute the matrices of elliptic projection on a face embedded in 3-D    
    Nv = size(P,1)
    ## derive local coordinates ne, te
    e1 = P[2, :] - P[1, :]
    en = P[1, :] - P[end, :]
    he = norm(e1)
    nf = cross(e1, en)
    nf ./= norm(nf)
    Ne = - cross(e1, nf)
    ne = Ne ./ fill(he, 3)  
    te = e1 ./ fill(he, 3)
    # transformation matrix
    Tmat = vcat(ne', te')
    # node for the 2-D polygon
    node = (P .- P[1, :]') / Tmat
    ## get auxiliary data
    verts = node
    verts1 = vcat(verts[2:end, :], verts[1, :]')
    area_components = verts[:, 1] .* verts1[:, 2] -
                      verts1[:,1] .* verts[:,2]
    ar = 0.5 * abs(sum(area_components))
    centroid = sum((verts + verts1) .* repeat(area_components, 1, 2), 
                   dims = 1) / (6*ar)
    diameter = maximum(norm(verts[l, :] - verts[m, :]) for l = 1:Nv, m = 1:Nv)
    ## compute the elliptic projection matrices
    # --- element information (only one element) ---
    xK, yK = centroid
    hK = diameter
    x, y = node[:, 1], node[:, 2]  
    # --- transition matrix ---
    D = hcat(1 .+ 0 * x, (x .- xK) / hK, (y .- yK) / hK)
    # --- elliptic projection ---
    # first term  = 0
    # second term
    rotid1 = [Nv, collect(1:(Nv - 1))...]
    rotid2 = [collect(2:Nv)..., 1] # ending and starting indices
    Gradm = [0 0; 1/hK 0; 0 1/hK]
    normVec = 0.5 * hcat(y[rotid2] - y[rotid1], x[rotid1] - x[rotid2])'
    # a rotation of edge vector
    B = Gradm * normVec
    # constraint
    Bs = B  
    Bs[1, :] .= 1 / Nv
    # consistency relation
    Gs = Bs * D
    # elliptic projection matrices
    Pis = Gs \ Bs  #Pi = D*Pis;
    return Pis
end
#-------------------------------------------------------------------------------
function faceTriangulation(elemf)
    # compute the triangulation of the faces of a polyhedra elemf
    # number of faces
    nf = length(elemf)
    # numbers of vertices on each face
    faceLen = length.(elemf) 
    # global index of vertices
    catfaces = Int64[]
    for f ∈ elemf
        append!(catfaces, f)
    end
    index3 = sort(unique(catfaces))
    totalid = [findfirst(isequal(catfaces[k]), index3) for 
                                                        k = 1:length(catfaces)]
    # local index of elemf
    elemfLocal, ck = Vector{Int64}[], 1
    for l ∈ faceLen
        push!(elemfLocal, totalid[ck:(ck + l - 1)])
        ck += l
    end
    # triangulation
    nTri = sum(faceLen .- 2)
    Tri = zeros(Int64, nTri, 3)  
    TriLocal = zeros(Int64, nTri, 3)
    id = 0
    for s = 1:nf
        face = elemf[s]
        faceLocal = elemfLocal[s]
        nt = faceLen[s] - 2
        Tri[(id + 1):(id + nt), :] .= hcat(fill(face[1], nt), face[2:(end-1)],
                                           face[3:end])
        TriLocal[(id + 1):(id + nt), :] .= hcat(fill(faceLocal[1], nt), 
                                                faceLocal[2:(end - 1)],
                                                faceLocal[3:end])
        id += nt
    end
    return Tri, index3, TriLocal, elemfLocal
end
#-------------------------------------------------------------------------------
function quadpts3(order)
    ## QUADPTS3 quadrature points in 3-D.
    #
    # [lambda,weight] = quadpts3(order) return quadrature points in a
    # tetrahedron with given order (up to 5) in the barycentric coordinates.
    #
    # The output lambda is a matrix of size nQ by 3, where nQ is the number of
    # quadrature points. lambda(i,:) is three barycentric coordinate of the
    # i-th quadrature point and lambda(:,j) is the j-th barycentric coordinate
    # of all quadrature points. The x-y-z coordinate of the p-th quadrature point
    # can be computed as 
    #
    #    pxyz = lambda(p,1)*node(elem(:,1),:) 
    #         + lambda(p,2)*node(elem(:,2),:)  
    #         + lambda(p,3)*node(elem(:,3),:)  
    #         + lambda(p,4)*node(elem(:,4),:);
    #
    # The weight of p-th quadrature point is given by weight(p). See
    # verifyquadpts for the usage of qudrature rules to compute integrals over
    # triangles.
    # 
    # References: 
    #
    # * Jinyun, Y. Symmetric Gaussian quadrature formulae for tetrahedronal
    # regions. Comput. Methods Appl. Mech. Engrg.. 43(3):349--353, 1984.
    #  
    # See also quadpts1, quadpts, verifyquadpts
    #
    # Copyright (C) Long Chen. See COPYRIGHT.txt for details.  
    if order > 5
        order = 5
    end
    lambda, weight = if order == 1 # order 1, nQuad 1
        lambda = [1/4, 1/4, 1/4, 1/4]
        weight = 1
        lambda, weight
    elseif order == 2 # order 2, nQuad 4
        alpha = 0.5854101966249685; 
        beta =  0.138196601125015;
        lambda = [alpha beta beta beta; 
                  beta alpha beta beta; 
                  beta beta alpha beta; 
                  beta beta beta alpha]
        weight = [1/4, 1/4, 1/4, 1/4]
        lambda, weight
    elseif order == 3 # order 3, nQuad 5
        lambda = [1/4 1/4 1/4 1/4; 
                  1/2 1/6 1/6 1/6; 
                  1/6 1/2 1/6 1/6; 
                  1/6 1/6 1/2 1/6; 
                  1/6 1/6 1/6 1/2]
        weight = [-4/5, 9/20, 9/20, 9/20, 9/20]
        lambda, weight
    elseif order == 4    # order 4, nQuad 16
        alpha1 = 0.7716429020672371
        beta1 =  0.7611903264425430e-1
        w1 = 0.5037379410012282e-1
        alpha = 0.4042339134672644
        beta = 0.7183164526766925e-1
        gamma = 0.11970052777978019
        w2 = 0.6654206863329239e-1
        lambda = [alpha1 beta1 beta1 beta1; 
                   beta1 alpha1 beta1 beta1; 
                   beta1 beta1 alpha1 beta1; 
                   beta1 beta1 beta1 alpha1; 
                   alpha alpha beta gamma; 
                   alpha alpha gamma beta; 
                   alpha beta alpha gamma; 
                   alpha beta gamma alpha; 
                   alpha gamma beta alpha; 
                   alpha gamma alpha beta; 
                   beta alpha alpha gamma; 
                   beta alpha gamma alpha; 
                   beta gamma alpha alpha; 
                   gamma alpha alpha beta; 
                   gamma alpha beta alpha; 
                   gamma beta alpha alpha]               
        weight = [w1, w1, w1, w1, w2, w2, w2, w2, w2, w2, w2, w2, w2, 
                  w2, w2, w2]
        lambda, weight
    elseif order == 5 # order 5, nQuad 17
        alpha1 = 0.7316369079576180
        beta1 =  0.8945436401412733e-1
        w1 = 0.6703858372604275e-1
        alpha = 0.4214394310662522
        beta = 0.2454003792903000e-1
        gamma = 0.1325810999384657
        w2 = 0.4528559236327399e-1
        lambda = [1/4 1/4 1/4 1/4; 
                  alpha1 beta1 beta1 beta1; 
                  beta1 alpha1 beta1 beta1; 
                  beta1 beta1 alpha1 beta1; 
                  beta1 beta1 beta1 alpha1; 
                  alpha alpha beta gamma; 
                  alpha alpha gamma beta; 
                  alpha beta alpha gamma; 
                  alpha beta gamma alpha; 
                  alpha gamma beta alpha; 
                  alpha gamma alpha beta; 
                  beta alpha alpha gamma; 
                  beta alpha gamma alpha; 
                  beta gamma alpha alpha; 
                  gamma alpha alpha beta; 
                  gamma alpha beta alpha; 
                  gamma beta alpha alpha]                  
        weight = [0.1884185567365411, w1, w1, w1, w1, w2, w2, w2, w2, w2, w2, 
                  w2, w2, w2, w2, w2, w2]
        lambda, weight
    end
    return lambda, weight
end    
#-------------------------------------------------------------------------------
function quadpts(order)
    ## QUADPTS quadrature points in 2-D.
    #
    # [lambda,weight] = quadpts(order) return quadrature points with given
    # order (up to 9) in the barycentric coordinates.
    #
    # The output lambda is a matrix of size nQ by 3, where nQ is the number of
    # quadrature points. lambda(i,:) is three barycentric coordinate of the
    # i-th quadrature point and lambda(:,j) is the j-th barycentric coordinate
    # of all quadrature points. The x-y coordinate of the p-th quadrature point
    # can be computed as 
    #
    #     pxy = lambda(p,1)*node(elem(:,1),:) 
    #         + lambda(p,2)*node(elem(:,2),:)  
    #         + lambda(p,3)*node(elem(:,3),:);
    #
    # The weight of p-th quadrature point is given by weight(p). See
    # verifyquadpts for the usage of qudrature rules to compute integrals over
    # triangles.
    # 
    # References: 
    #
    # * David Dunavant. High degree efficient symmetrical Gaussian
    #    quadrature rules for the triangle. International journal for numerical
    #    methods in engineering. 21(6):1129--1148, 1985. 
    # * John Burkardt. DUNAVANT Quadrature Rules for the Triangle.
    #    http://people.sc.fsu.edu/~burkardt/m_src/dunavant/dunavant.html
    # 
    # See also quadpts1, quadpts3, verifyquadpts
    #
    # Order 6 - 9 is added by Huayi Wei, modify by Jie Zhou
    #
    # Copyright (C) Long Chen. See COPYRIGHT.txt for details. 
    if order > 9
        order = 9
    end
    lambda, weight = if order == 1     # order 1, nQuad 1
        lambda = [1/3, 1/3, 1/3]
        weight = 1
        lambda, weight
    elseif order == 2 # order 2, nQuad 3
        lambda = [2/3 1/6 1/6; 
                  1/6 2/3 1/6; 
                  1/6 1/6 2/3]
        weight = [1/3 1/3 1/3]
        lambda, weight
    elseif  order == 3 # order 3, nQuad 4
        lambda = [1/3 1/3 1/3; 
                  0.6 0.2 0.2; 
                  0.2 0.6 0.2; 
                  0.2 0.2 0.6]
        weight = [-27/48, 25/48, 25/48, 25/48]
        lambda, weight
    elseif order == 4 # order 4, nQuad 6
        lambda = [0.108103018168070 0.445948490915965 0.445948490915965; 
                  0.445948490915965 0.108103018168070 0.445948490915965; 
                  0.445948490915965 0.445948490915965 0.108103018168070; 
                  0.816847572980459 0.091576213509771 0.091576213509771; 
                  0.091576213509771 0.816847572980459 0.091576213509771; 
                  0.091576213509771 0.091576213509771 0.816847572980459]
        weight = [0.223381589678011, 0.223381589678011, 0.223381589678011, 
                  0.109951743655322, 0.109951743655322, 0.109951743655322]
        lambda, weight
    elseif order == 5 # order 5, nQuad 7
        alpha1 = 0.059715871789770;      beta1 = 0.470142064105115
        alpha2 = 0.797426985353087;      beta2 = 0.101286507323456
        lambda = [   1/3    1/3    1/3; 
                   alpha1  beta1  beta1; 
                    beta1 alpha1  beta1; 
                    beta1  beta1 alpha1; 
                   alpha2  beta2  beta2; 
                    beta2 alpha2  beta2; 
                    beta2  beta2 alpha2]
        weight = [0.225, 0.132394152788506, 0.132394152788506, 
                  0.132394152788506, 0.125939180544827, 0.125939180544827, 
                  0.125939180544827]
        lambda, weight
    elseif order == 6        
        A =[0.249286745170910  0.249286745170910  0.116786275726379
            0.249286745170910  0.501426509658179  0.116786275726379
            0.501426509658179  0.249286745170910  0.116786275726379
            0.063089014491502  0.063089014491502  0.050844906370207
            0.063089014491502  0.873821971016996  0.050844906370207
            0.873821971016996  0.063089014491502  0.050844906370207
            0.310352451033784  0.636502499121399  0.082851075618374
            0.636502499121399  0.053145049844817  0.082851075618374
            0.053145049844817  0.310352451033784  0.082851075618374
            0.636502499121399  0.310352451033784  0.082851075618374
            0.310352451033784  0.053145049844817  0.082851075618374
            0.053145049844817  0.636502499121399  0.082851075618374]
        lambda = hcat(A[:, [1, 2]], 1 .- sum(A[:,[1, 2]], dims = 2))
        weight = A[:, 3]
        lambda, weight
    elseif order == 7
        A =[0.333333333333333  0.333333333333333 -0.149570044467682
            0.260345966079040  0.260345966079040  0.175615257433208
            0.260345966079040  0.479308067841920  0.175615257433208
            0.479308067841920  0.260345966079040  0.175615257433208
            0.065130102902216  0.065130102902216  0.053347235608838
            0.065130102902216  0.869739794195568  0.053347235608838
            0.869739794195568  0.065130102902216  0.053347235608838
            0.312865496004874  0.638444188569810  0.077113760890257
            0.638444188569810  0.048690315425316  0.077113760890257
            0.048690315425316  0.312865496004874  0.077113760890257
            0.638444188569810  0.312865496004874  0.077113760890257
            0.312865496004874  0.048690315425316  0.077113760890257
            0.048690315425316  0.638444188569810  0.077113760890257]
        lambda = hcat(A[:, [1, 2]], 1 .- sum(A[:,[1, 2]], dims = 2))
        weight = A[:, 3]
        lambda, weight
    elseif order == 8
        A =[0.333333333333333  0.333333333333333  0.144315607677787
            0.081414823414554  0.459292588292723  0.095091634267285
            0.459292588292723  0.081414823414554  0.095091634267285
            0.459292588292723  0.459292588292723  0.095091634267285
            0.658861384496480  0.170569307751760  0.103217370534718
            0.170569307751760  0.658861384496480  0.103217370534718
            0.170569307751760  0.170569307751760  0.103217370534718
            0.898905543365938  0.050547228317031  0.032458497623198
            0.050547228317031  0.898905543365938  0.032458497623198
            0.050547228317031  0.050547228317031  0.032458497623198
            0.008394777409958  0.263112829634638  0.027230314174435
            0.008394777409958  0.728492392955404  0.027230314174435
            0.263112829634638  0.008394777409958  0.027230314174435
            0.728492392955404  0.008394777409958  0.027230314174435
            0.263112829634638  0.728492392955404  0.027230314174435
            0.728492392955404  0.263112829634638  0.027230314174435]
        lambda = hcat(A[:, [1, 2]], 1 .- sum(A[:,[1, 2]], dims = 2))
        weight = A[:, 3]     
        lambda, weight   
    elseif order == 9
        A =[0.333333333333333  0.333333333333333  0.097135796282799
            0.020634961602525  0.489682519198738  0.031334700227139
            0.489682519198738  0.020634961602525  0.031334700227139
            0.489682519198738  0.489682519198738  0.031334700227139
            0.125820817014127  0.437089591492937  0.07782754100474
            0.437089591492937  0.125820817014127  0.07782754100474
            0.437089591492937  0.437089591492937  0.07782754100474
            0.623592928761935  0.188203535619033  0.079647738927210
            0.188203535619033  0.623592928761935  0.079647738927210
            0.188203535619033  0.188203535619033  0.079647738927210
            0.910540973211095  0.044729513394453  0.025577675658698
            0.044729513394453  0.910540973211095  0.025577675658698
            0.044729513394453  0.044729513394453  0.025577675658698
            0.036838412054736  0.221962989160766  0.043283539377289
            0.036838412054736  0.741198598784498  0.043283539377289
            0.221962989160766  0.036838412054736  0.043283539377289
            0.741198598784498  0.036838412054736  0.043283539377289
            0.221962989160766  0.741198598784498  0.043283539377289
            0.741198598784498  0.221962989160766  0.043283539377289]
        lambda = hcat(A[:, [1, 2]], 1 .- sum(A[:,[1, 2]], dims = 2))
        weight = A[:, 3]
        lambda, weight
    end
    return lambda, weight
end
#-------------------------------------------------------------------------------
function integralPolyhedron(fun, n, node3, elemf, centroid3)
    # Approximate integrals in a polyhedral domain
    # n: n-th order quadrature rule
    # fun: one or more anonymous functions, e.g. 
    #    fun = @(x,y,z) [f1(x,y,z), f2(x,y,z)]
    ## Triangulation and geo of the polyhedon
    # triangulation of faces
    T = typeof(node3[1])
    Tri, index3, TriLocal = PolygonalFem.faceTriangulation(elemf)
    V = node3[index3, :]
    if size(Tri, 1) > 4
        # triangulation of the polyhedron    
        nodeTet = vcat(V, centroid3')
        elemTet = hcat(TriLocal, fill(maximum(TriLocal) + 1, size(TriLocal,1)))
    else # tetrahedron
        nodeTet = V
        elemTet = collect(1:4)'
    end
    # volumes of all tetrahedrons
    volume, _ = simplexvolume(nodeTet, elemTet)
    ## Guass-Quadrature
    lambda, weight = PolygonalFem.quadpts3(n)
    nQuad = length(weight)
    sout = length(fun(rand(3)))
    Int = [zeros(T, sout) for _ = 1:size(elemTet, 1)]
    for p = 1:nQuad
        pz = lambda[p, 1] * nodeTet[elemTet[:, 1], :] +
             lambda[p, 2] * nodeTet[elemTet[:, 2], :] +
             lambda[p, 3] * nodeTet[elemTet[:, 3], :] +
             lambda[p, 4] * nodeTet[elemTet[:, 4], :]
        fp = [fun(pz[k, :]) for k = 1:size(pz, 1)]
        Int .+= weight[p] * fp
    end
   return sum(volume .* Int)
end
function simplexvolume(node, elem)
    NT = size(elem, 1)
    # tetrahedra
    d12 = node[elem[:, 2], :] - node[elem[:, 1], :]
    d13 = node[elem[:, 3], :] - node[elem[:, 1], :]
    d14 = node[elem[:, 4], :] - node[elem[:, 1], :]
    v = [dot(cross(d12[k, :], d13[k, :]), d14[k, :]) / 6 for k = 1:NT]
    idx = findall(v .< 0) 
    v[idx, :] = - v[idx, :]
    elemSign = ones(Int64, NT)
    elemSign[idx] .= -1
    return v, elemSign 
end
#-------------------------------------------------------------------------------
function polyarea(P)
    # Compute the area of a polygon in 3-D
    nv = size(P, 1)
    tri = hcat(ones(Int64, nv - 2), collect(2:(nv - 1)), collect(3:nv))
    d12 = P[tri[:, 2], :] - P[tri[:, 1], :]
    d13 = P[tri[:, 3], :] - P[tri[:, 1], :]
    nbtri = size(tri, 1)
    normal = [cross(d12[k, :], d13[k, :]) for k ∈ 1:nbtri]
    areaTri = [norm(normal[k]) for k ∈ 1:nbtri] / 2
    area = sum(areaTri)
    return area
end
########################### BUGGED NOT USED ####################################
function buggy_polycentroid3(V, Tri)
    # polycentroid3 returns the x,y,z coordinates of centroid
    # of surface triangulated polyhedron.
    #
    # INPUT:
    #     vertex: Point Cloud of Shape
    #         V(:,1) : x coordinates
    #         V(:,2) : y coordinates
    #         V(:,3) : z coordinates
    #     Tri: connectivity list of face triangulation
    vector1 = V[Tri[:, 2], :] - V[Tri[:, 1], :]
    vector2 = V[Tri[:, 3], :] - V[Tri[:, 1], :]
    nTri = size(Tri, 1)
    triangAreasTmp = 0.5 * [cross(vector1[k, :], vector2[k, :]) for k = 1:nTri]
    triangAreas = norm.(triangAreasTmp) # area of each triangle
    totArea = sum(triangAreas) # total area
    point1 = V[Tri[:, 1], :]
    point2 = V[Tri[:, 2], :]
    point3 = V[Tri[:, 3], :]
    centroidTri = 1/3 * (point1 + point2 + point3) # cent. of each triangle
    mg2 = triangAreas .* centroidTri
    centroid = sum(mg2, dims = 1)[:] ./ totArea
    return centroid
end