export
AABB,
overlaps,
contains,
inflate,
volume,
refit,
randAABB,
centre,
radius,
radius_squared

using Base.Iterators
using LinearAlgebra
import Base.union
import Base.isless

struct AABB{T}
  min::Vector3{T}
  max::Vector3{T}
end

AABB{T}() where T = AABB{T}(ones(Vector3{T}), ones(Vector3{T}))

function AABB( minx::T, maxx::T, miny::T, maxy::T, minz::T, maxz::T ) where {T}
  AABB(Vector3([minx, miny, minz]), Vector3([maxx, maxy, maxz]))
end

AABB( P::Vector{Point{T}} ) where {T} = reduce(union,(p->AABB(p,p)).(P)) 

function AABB( a::Point{T}, b::Point{T}, c::Point{T} ) where {T}
  @inbounds r = AABB(min(a[1], b[1], c[1]),
                     max(a[1], b[1], c[1]),
                     min(a[2], b[2], c[2]),
                     max(a[2], b[2], c[2]),
                     min(a[3], b[3], c[3]),
                     max(a[3], b[3], c[3]))
  r
end

function AABB( a::Vector3{T}, b::Vector3{T}, c::Vector3{T}, rest::Vector3{T}...) where {T} 
  @inbounds r = AABB(min(a[1],b[1],c[1],map(p->p[1],rest)...), max(a[1],b[1],c[1],map(p->p[1],rest)...), 
                     min(a[2],b[2],c[2],map(p->p[2],rest)...), max(a[2],b[2],c[2],map(p->p[2],rest)...), 
                     min(a[3],b[3],c[3],map(p->p[3],rest)...), max(a[3],b[3],c[3],map(p->p[3],rest)...)); r
end

function refit( aabb::AABB{T} ) where T
  a = aabb.min
  b = aabb.max
  @inbounds AABB(Vector3(min(a[1],b[1]), min(a[2],b[2]), min(a[3],b[3])),
                 Vector3(max(a[1],b[1]), max(a[2],b[2]), max(a[3],b[3])))
end

randAABB(r = 0:0.001:1) = refit(AABB(rand(r), rand(r), rand(r), rand(r), rand(r), rand(r)))
Base.rand(::Type{AABB{Float64}}) = randAABB()

"""
"""
function union( a::AABB{T}, b::AABB{T} ) where {T}
  @inbounds @fastmath r = AABB(min(a.min[1],b.min[1]), max(a.max[1],b.max[1]),
                               min(a.min[2],b.min[2]), max(a.max[2],b.max[2]),
                               min(a.min[3],b.min[3]), max(a.max[3],b.max[3])); r
end

"""
check if all extents of a are < all extents of b.  if so, then disjoint.
"""
function isless(a::AABB{T}, b::AABB{T}) where {T}
  @inbounds @fastmath r = (a.max[1] < b.min[1]) || (a.max[2] < b.min[2]) || (a.max[3] < b.min[3]); r
end


function overlaps(a::AABB{T}, b::AABB{T}) where {T}
  @inbounds @fastmath r = !((a < b) || (b < a)); r
end


function contains(a::AABB{T}, b::AABB{T}) where {T}
  @inbounds @fastmath r = ((a.max[1] >= b.max[1]) && (a.min[1] <= b.min[1])) &&
  ((a.max[2] >= b.max[2]) && (a.min[2] <= b.min[2])) &&
  ((a.max[3] >= b.max[3]) && (a.min[3] <= b.min[3])); r
end


function inflate(a::AABB{T}, r::S) where {T,S}
  @inbounds @fastmath rv = Vector3([r,r,r])
  @inbounds @fastmath r = AABB( a.min-rv, a.max+rv ); r
end


function volume(a::AABB{T}) where {T}
  @inbounds @fastmath r = (a.max[1]-a.min[1])*(a.max[2]-a.min[2])*(a.max[3]-a.min[3]); r
end

function centre(a::AABB{T}) where {T}
  @inbounds @fastmath (a.max+a.min)*T(0.5)
end

function radius(a::AABB{T}) where {T}
  @inbounds @fastmath norm(a.max-a.min)*T(0.5)
end

function radius_squared(a::AABB{T}) where {T}
  @inbounds @fastmath LinearAlgebra.norm_sqr(a.max-a.min)*T(0.5)
end

