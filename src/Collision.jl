"""
Collision detection
"""
module Collision

export triangle_triangle, triangle_line, triangle_plane, plane_segement, γpoint

using StaticArrays, SparseArrays, LinearAlgebra, Base.Iterators, IterTools

partial = (f::Function,y...)->(z...)->f(y...,z...)
second(x) = x[2]

# create a basis for Clifford Algebra using Gamma Matrices
# we will do some projective geometry with it.
# ∧ product works pretty well, but regressive product doesn't.  That's ok we just need join
const 𝑖 = 1.0im
γ0 = [1.0 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 -1]
γ1 = [0.0 0 0 1; 0 0 1 0; 0 -1 0 0; -1 0 0 0]
γ2 = [0.0 0 0 -𝑖; 0 0 𝑖 0; 0 𝑖 0 0; -𝑖 0 0 0]
γ3 = [0.0 0 1 0; 0 0 0 -1; -1 0 0 0; 0 1 0 0]

# (-+++) clifford algebra
# should do a performance comparison with StaticArrays, SparseArrays and Arrays
const g1,g2,g3,g4 = sparse.(im .* (γ0, γ1, γ2, γ3))

const 𝐼 = g1*g2*g3*g4
const ProjectivePoint = typeof(g1)
const ProjectiveLine = ProjectivePoint
const ProjectivePlane = ProjectivePoint

∧(a::T, b::T) where T<:SparseMatrixCSC = a*b
const join = ∧

# orientation of a pseudoscalar element.
import Base.sign
sign(ps) = sign(real(tr(ps*𝐼)))

""" 
    γpoint(x,y,z) 

create a homogeneous point for collision detection using our specialized basis.
simply pass in the euclidean x,y,z coordinates
"""
γpoint(x,y,z) = x*g2 + y*g3 + z*g4 + 1.0g1
γpoint(xyz) = γpoint(xyz...)

""" true if triangle and line intersect """
function triangle_line((a,b,c)::Tuple{P,P,P}, (d,e)::Tuple{P,P}) where P <: ProjectivePoint
  de = d∧e
  crossings = map(sign, (de∧a∧b, de∧b∧c, de∧c∧a))

  # line passes through inside of triangle if line-crossing volumes are all same orientation
  # if product is zero then we are coplanar or exactly on an edge 
  sum(crossings)^2 == 9.0 || prod(crossings) == 0.0
end

""" true if triangle and plane intersect """
triangle_plane((a,b,c)::Tuple{P,P,P}, p::V) where {V<:ProjectivePlane, P<:ProjectivePoint} =
  !(sign(p∧a) == sign(p∧b) == sign(p∧c))


""" true if triangle and segment intersect """
plane_segment(abc::V, (d,e)::Tuple{P,P}) where {V<:ProjectivePlane, P<:ProjectivePoint} =
  sign(abc∧d) != sign(abc∧e)

""" 
    triangle_triangle((a,b,c), (d,e,f))

return true if the triangle defined by points a,b,c is intersecting the triangle defined by points d,e,f
operates on projective points created by γpoint
"""
function triangle_triangle((a,b,c)::T, (d,e,f)::T) where {P<:ProjectivePoint, T<:Tuple{P,P,P}}
  abc = a∧b∧c
  def = d∧e∧f
  abc_edges = ((a,b), (b,c), (c,a))
  def_edges = ((d,e), (e,f), (f,d))

  check_plane(plane, edges) = ( plane_segment(plane, edge) for edge in edges )
  check_plane2(plane, edges) = imap(second, Iterators.filter(first, zip(check_plane(plane, edges), edges)))

  check_tri(tri, plane, edges) = ( triangle_line(tri, edge) for edge in check_plane2(plane, edges) )

  for ihit in check_tri((a,b,c), abc, def_edges)
    if ihit
      return true
    end
  end
  for ihit in check_tri((d,e,f), def, abc_edges)
    if ihit
      return true
    end
  end

  return false
end

function triangle_triangle((a,b,c)::T, (d,e,f)::T) where {P<:SVector, T<:Tuple{P,P,P}}
  triangle_triangle((γpoint(a), γpoint(b), γpoint(c)), (γpoint(d), γpoint(e), γpoint(f)))
end

triangle_triangle(a,b,c,d,e,f) = triangle_triangle((a,b,c), (d,e,f))

end
