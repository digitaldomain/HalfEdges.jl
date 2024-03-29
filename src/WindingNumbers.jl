det22(a,b,c,d) = a*d-c*b
det33(a,b,c) = @inbounds a[1]*det22(b[2],c[2],b[3],c[3]) - b[1]*det22(a[2],c[2],a[3],c[3]) + c[1]*det22(a[2],b[2],a[3],b[3])

"""
  solid_angle(topo, P, he, p)

  
the solid angle for the signed surface area of the triangle of he projected onto unit sphere centred at p 
"""
function solid_angle(topo::Topology, P::TA, h::HalfEdgeHandle, p::T) where {T, TA<:Vector{T}}
  if isboundary(topo, h) 
    return 0.0
  end

  @inbounds a = P[head(topo, h)] - p
  @inbounds b = P[head(topo, next(topo, h))] - p
  @inbounds c = P[head(topo, next(topo, next(topo, h)))] - p
  aa = norm(a); bb = norm(b); cc = norm(c)
  2.0*atan(det33(a,b,c), aa*bb*cc + a⋅b*cc + b⋅c*aa + c⋅a*bb)
end

function solid_angle_normal(topo::Topology, P::TA, h::HalfEdgeHandle, p::T) where {T, TA<:Vector{T}}
  if isboundary(topo, h) 
    return 0.0
  end

  @inbounds a = P[head(topo, h)] - p
  @inbounds b = P[head(topo, next(topo, h))] - p
  @inbounds c = P[head(topo, next(topo, next(topo, h)))] - p
  aa = norm(a); bb = norm(b); cc = norm(c)
  Ω = 2.0*atan(det33(a,b,c), aa*bb*cc + a⋅b*cc + b⋅c*aa + c⋅a*bb)


  if isapprox(Ω, 0.0, atol=1e-16)
    zero(T) 
  else
    w = LinearAlgebra.norm(a+b+c)/3.0
    centroidn = normalize((normalize(a) + normalize(b) + normalize(c)))
    centroidn*Ω*w
  end
end

function winding_normal(topo::Topology, P, q::VT) where {T<:Number, 
                                                         VT<:AbstractVector{T}}
  mapreduce(i->solid_angle_normal(topo, P, i, q), +, faces(topo))/(4.0π)
end


#=== some tensor stuff ===#
# need this to make approximate winding numbers managable
# works pretty seamlessly with julia's multidimensional arrays. 
# pops out an Array quite easily from pretty much any sensible operation
# make sure to convert to Array after precomputation as this isn't very performant for access

struct LazyTensor{T,N} <: AbstractArray{T,N}
  V::Union{NTuple{N}, Vector{LazyTensor{T,N}}}
  dim
end

Base.size(A::LazyTensor) = A.dim
Base.getindex(A::LazyTensor{T,N}, I::Vararg{Int,NN}) where {T,N,NN} = getlazyindex(A.V, I...)

getlazyindex(V::NTuple{N}, I::Vararg{Int,NN}) where {N,NN} = prod(map(i->V[i][I[i]], 1:length(I)))

function getlazyindex(V::Vector{LazyTensor{T,N}}, I::Vararg{Int,NN}) where {T,N,NN} 
  sum(map(Aᵢ->getindex(Aᵢ, I...), V))
end

Base.IndexStyle(::Type{LazyTensor}) = IndexCartesian()

cons(a::T, b::NTuple{N,T}) where {N,T} = (a, b...) 
cons(a::T, b::Vector{T}) where T = vcat(a, b)
cons(a::NTuple{N,T}, b::T) where {N,T} = (a..., b) 
cons(a::Vector{T}, b::T) where T = vcat(a, b)

""" distributes ⊗ into a + b +... """
cons(a::Vector{LazyTensor{T,N}}, b::AbstractVector{T}) where {N,T} = map(aᵢ->aᵢ⊗b, a)
cons(a::AbstractVector{T}, b::Vector{LazyTensor{T,N}}) where {N,T} = map(bᵢ->a⊗bᵢ, b)

⊗(a::T, b::T) where {F<:Real, T<:AbstractVector{F}} = LazyTensor{F,2}((a,b), (length(a), length(b)))

⊗(a::T, b::VT) where {F<:Real, T<:AbstractVector{F}, N, VT<:LazyTensor{F,N}} = 
  LazyTensor{F,N+1}(cons(a, b.V), (length(a), b.dim...))

⊗(a::VT, b::T) where {F<:Real, T<:AbstractVector{F}, N, VT<:LazyTensor{F,N}} = 
  LazyTensor{F,N+1}(cons(a.V, b), (a.dim..., length(b)))

#!me don't actually need this, so didn't finish implementing
#⊗(a::TA, b::TB) where {F<:Real, NA, TA<:LazyTensor{F,NA}, NB, TB<:LazyTensor{F,NB}} = 
#  LazyTensor{F,NA+NB}(cons(a.V, b.V), (a.dim..., b.dim...))

⊗(a) = ⊗(a, a)

Base.:+(a::LT, b::LT) where {T, N, LT<:LazyTensor{T,N}} = LazyTensor{T, N}([a,b], a.dim)


# dipole derivatives from appendix A of "Fast Winding Numbers for Soups and Clouds"
function ∇G(q,x) 
  r = x-q
  r/(4.0π*norm(r)^3)
end

function ∇²G(q,x)
  r = x-q
  rlen = norm(r)
  I/(4.0π*rlen^3) - 3.0*(r*r')/(4.0π*rlen^5)
end

# literally from the paper, lots of 0*0 happening and that ⊗ is brutal slow
# likely we wont actually use this third order expansion, so don't worry about optz right now
function ∇³G(q,x)
  r = x-q
  rlen = norm(r)
  Array(sum(map(eᵢ->r⊗eᵢ⊗eᵢ + eᵢ⊗r⊗eᵢ + eᵢ⊗eᵢ⊗r, [SVector(1.0,0.0,0.0),
                                                  SVector(0.0,1.0,0.0),
                                                  SVector(0.0,0.0,1.0)])))/(-4.0π*rlen^5) + 
  15.0*Array(r⊗r⊗r)/(4.0π*rlen^7)
end

function centroid(topo, P,  F::Vector{HalfEdgeHandle})
  Pees = unique(mapreduce(partial(vertices, topo), vcat, F; init=[]))
  p̃ = sum(P[Pees])/length(Pees)
end

# triangle integrals from appendix B of "Fast Winding Numbers for Soups and Clouds"

struct ApproxWindingNumber2
  p̃::SVector{3, Float64}
  term1::SVector{3, Float64}
  term2::SMatrix{3,3, Float64}
end

struct ApproxWindingNumber3
  p̃::SVector{3, Float64}
  term1::SVector{3, Float64}
  term2::SMatrix{3,3, Float64}
  term3::Array{Float64, 3}
end

approx_winding_number(awn::ApproxWindingNumber2, 
                      q::SVector{3,Float64}) = awn.term1⋅∇G(q, awn.p̃) + awn.term2⋅∇²G(q, awn.p̃)

"""
    approx_winding_number(topo, P, F, order=3)

assumes triangle mesh.  approximate winding number function for a cluster of triangles
"""
function approx_winding_number(topo, P, F::Vector{HalfEdgeHandle}, order=2)

  p̃ = centroid(topo, P, F)

  #!me check https://github.com/alecjacobson/WindingNumber/blob/master/UT_SolidAngle.cpp 
  #!me says we need to use NORMALIZED normal to multiply integrals by, paper says otherwise...
  ant = map(partial(weightednormal, topo, P), F)

  term1 = sum(ant)
  term2 = map(zip(ant, F)) do (antᵢ, heh)
    (sum(P[vertices(topo, heh)])/3.0 - p̃)⊗antᵢ
  end |> sum |> SMatrix{3,3}∘Array

  # very inefficient, but order 2 seems good enough, so don't worry about it for now
  if order == 3
    term3 = map(zip(ant, F)) do (antᵢ, heh)
      xᵢ, xⱼ, xk = P[vertices(topo, heh)]
      Ct = ⊗(0.5*(xᵢ+xⱼ)-p̃) +
           ⊗(0.5*(xⱼ+xk)-p̃) +
           ⊗(0.5*(xk+xᵢ)-p̃)
      Ct⊗antᵢ
    end |> sum |> partial(*, 1.0/6.0) |> Array

  #  w̃ = (q)->term1⋅∇G(q, p̃) + term2⋅∇²G(q, p̃) + term3⋅∇³G(q, p̃)
    ApproxWindingNumber3(p̃, term1, term2, term3) 
  else
  #  w̃ = (q)->term1⋅∇G(q, p̃) + term2⋅∇²G(q, p̃)
    ApproxWindingNumber2(p̃, term1, term2) 
  end

  # is storing a closure reasonably performant?  if not we should just be storing (p̃, term1...)
  #return w̃

end

approx_winding_number(topo, P, F::Vector{FaceHandle}) = approx_winding_number(topo, P, faces(topo)[F])

function approx_winding_number(topo, P, q, node::BT) where {BT<:BinaryTree}
  F = ileaves(node) |> collect |> c->map(value∘data,c)
  approx_winding_number(topo, P, map(i->halfedge(topo, FaceHandle(i)), F))(q)
end

struct WindingNumberCache
  bvh::IndexedBinaryTree{AVLData{AABBNodeData{Float64,Int64}}}
  #W::Vector{Function}
  W::Vector{ApproxWindingNumber2}  #!me ignore order 3 for now
  height_cutoff::Int64
  β::Float64
  F::Vector{Vector{HalfEdgeHandle}}
  WindingNumberCache(a,b,c,d,e) = new(a,b,c,d,e)
  WindingNumberCache(a,b,c,d) = new(a,b,c,d,Vector{Vector{HalfEdgeHandle}})
end

WindingNumberCache(a,b,c) = WindingNumberCache(a, b, c, 2.0)

function winding_number_cache(topo, P, tris_for_brute = 100, order=2)
  aabbs = map((Fᵢ, i)->AABBNodeData(AABB(P[Fᵢ]), i), polygons(topo), 1:nfaces(topo))
  bt = avlitree(map( (Fᵢ, i)->AABBNodeData(AABB(P[Fᵢ]), i), polygons(topo), 1:nfaces(topo)))

  # threshold for when we generate approx w
  height_for_approx = (minheight(tris_for_brute) + maxheight(tris_for_brute))/2 |> Int∘round

  #nnodes = mapreduce(index, max, TraverseWhen(n->height(n) >= height_for_approx, bt, identity))
  nnodes = mapreduce(index, max, itraversewhen(n->height(n) >= height_for_approx, bt, eltype(bt), identity))
  F = polygons(topo)

#  lc = TraverseWhen(n->height(n) >= height_for_approx, bt, Leaves)
  lc = TraverseWhen(n->height(n) >= height_for_approx, bt, ileaves)
  
  faceheh(n) = n |> data |> value |> i->halfedge(topo, FaceHandle(i))
  W = ApproxWindingNumber2[]
  resize!(W, nnodes)
  W[ TraverseWhen(n->height(n) >= height_for_approx, bt, index) |> collect ] = 
    [ approx_winding_number(topo, P, map(faceheh, lcᵢ), order) for lcᵢ in lc ]

  #  W[ itraversewhen(n->height(n) >= height_for_approx, bt, Int, index) ] = 
  #    [ approx_winding_number(topo, P, map(faceheh, lcᵢ), order) for lcᵢ in lc ]


  at_height = Iterators.filter(n->height(n) == height_for_approx-1,
                               itraversewhen(n->height(n) >= height_for_approx-1, 
                                             bt, eltype(bt), identity))
  iF = map(index, at_height)
  nF = reduce(max, iF) 
  F = Vector{HalfEdgeHandle}[]
  resize!(F, nF)
  F[iF] = map(n->faceheh.(ileaves(n)), at_height)

  WindingNumberCache(bt, W, height_for_approx, 2.0, F)
end

"""
"""
function winding_number(topo, P, q, cache::WindingNumberCache)
  node = cache.bvh
  β² = cache.β*cache.β
  winding_number(topo, P, q, node, β², cache)
end

function winding_number(topo, P, q, node::BT, β²::Float64, cache) where {T, BT<:IndexedBinaryTree{T}}
  if height(node) == cache.height_cutoff-1
    return mapreduce(heh->solid_angle(topo, P, heh, q), +, cache.F[index(node)])/(4.0π)
  elseif height(node) < cache.height_cutoff
    return winding_number(topo, P, q, node)
  else
    aabb = key(node)
    p = centre(aabb)
    r = radius_squared(aabb)
    if LinearAlgebra.norm_sqr(q - p) > β²*r
      #println("WHATUP", index(node))
      #return cache.W[index(node)](q)
      return approx_winding_number(cache.W[index(node)], q)
    else
      # descend 
      #println("crazy!!!")
      return winding_number(topo, P, q, left(node), β², cache) +
             winding_number(topo, P, q, right(node), β², cache)
    end
  end
end

function iwinding_number(topo, P, q, cache::WindingNumberCache)
  node = cache.bvh
  β² = cache.β*cache.β
  iwinding_number(topo, P, q, node, β², cache)
end

function iwinding_number(topo, P, q, node::BT, β²::Float64, cache) where {T, BT<:IndexedBinaryTree{T}}
  sum = 0.0
  nodestack = [index(node)] 

  while !isempty(nodestack)
    node = cache.bvh[pop!(nodestack)]
    #println(height(node), " compare ", cache.height_cutoff-1)
    if height(node) < cache.height_cutoff-1
      sum += winding_number(topo, P, q, node)
    elseif height(node) == cache.height_cutoff-1
      sum += mapreduce(heh->solid_angle(topo, P, heh, q), +, cache.F[index(node)])/(4.0π)
      #  sum += solid_angle(topo, P, halfedge(topo, FaceHandle(value(data(node)))), q)/(4.0π)
      #if height(node) == 0
      #  sum += solid_angle(topo, P, halfedge(topo, FaceHandle(value(data(node)))), q)/(4.0π)
    else
      aabb = key(node)
      p = centre(aabb)
      r = radius_squared(aabb)
      if LinearAlgebra.norm_sqr(q - p) > β²*r
        #println("WHATUP", index(node))
        #return cache.W[index(node)](q)
        sum += approx_winding_number(cache.W[index(node)], q)
      else
        # descend 
        #println("crazy!!!")
        push!(nodestack, ileft(node))
        push!(nodestack, iright(node))
      end
    end
  end
  sum
end

function winding_number(topo, P, q, node::BT) where {T, BT<:IndexedBinaryTree{T}}
  #mapreduce(l->solid_angle(topo, P, halfedge(topo, FaceHandle(value(data(l)))), q), +, ileaves(node))/(4.0π)
  wn = 0.0
  nodestack=[node]
  while !isempty(nodestack)
    node = pop!(nodestack)

    if isleaf(node)
      wn += solid_angle(topo, P, halfedge(topo, FaceHandle(value(data(node)))), q)
    else
      l = left(node)
      r = right(node)
      l != BoundingVolumeTrees.EmptyIndexedNode && push!(nodestack, left(node))
      r != BoundingVolumeTrees.EmptyIndexedNode && push!(nodestack, right(node))
    end
  end
  wn/(4.0π)
end

winding_number(topo, P, q) = mapreduce(i->solid_angle(topo, P, i, q), +, faces(topo))/(4.0π)

winding_numbers(topo, P::V, cache; jiggle = eps(T)) where {T, 
                                                           PT<:AbstractVector{T}, 
                                                           V<:AbstractVector{PT}} = 
  map( q->iwinding_number(topo, P, q, cache), P+vertexnormals(topo, P)*jiggle)

"""
    winding_numbers(topo, P::V; jiggle = ϵ)

The winding number for each vertex of a mesh.  
A vertex is lifted slightly off the surface in the normal direction to give winding number of 0 for vertices on the surface, but otherwise on exterior of the mesh.
The amount of normal lifting is controlled with the jiggle parameter 
"""
function winding_numbers(topo, P::V; jiggle = eps(T), approx=false, cacheleaves=50 ) where {T, PT<:AbstractVector{T}, V<:AbstractVector{PT}}
  if approx
    #!me fast winding numbers are anything but atm. some serious performance issue we need to investigate
    cache = winding_number_cache(topo, P, cacheleaves)
    winding_numbers(topo, P, cache; jiggle=jiggle)
    #map( q->iwinding_number(topo, P, q, cache), P+vertexnormals(topo, P)*jiggle)
  else
    map( q->winding_number(topo, P, q), P+vertexnormals(topo, P)*jiggle)
  end
end

