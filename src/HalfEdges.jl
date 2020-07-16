# Copyright 2020 Digital Domain 3.0
#
# Licensed under the Apache License, Version 2.0 (the "Apache License")
# with the following modification; you may not use this file except in
# compliance with the Apache License and the following modification to it:
# Section 6. Trademarks. is deleted and replaced with:
#
# 6. Trademarks. This License does not grant permission to use the trade
#    names, trademarks, service marks, or product names of the Licensor
#    and its affiliates, except as required to comply with Section 4(c) of
#    the License and to reproduce the content of the NOTICE file.
#
# You may obtain a copy of the Apache License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the Apache License with the above modification is
# distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied. See the Apache License for the specific
# language governing permissions and limitations under the Apache License.

module HalfEdges

using Base.Iterators, LinearAlgebra, SparseArrays

export 
  Topology,
  HalfEdge,
  HalfEdgeHandle,
  VertexHandle,
  FaceHandle,
  EdgeHandle,
  OneRing,
  Rim,
  Polygon,
  BoundaryLoop,
  UniqueHalfEdges,
  IncidentEdges,
  IncidentHalfEdges,
  nvertices,
  nverts,
  nedges,
  nfaces,
  nhalfedges,
  facelist,
  incidence,
  loadmesh,
  head,
  tail,
  isboundary,
  opposite,
  boundary_verts,
  boundary_vertices,
  boundary_interior,
  vertexnormals,
  normals,
  polygons,
  face,
  halfedge,
  vertices,
  winding_number,
  winding_number_cache

include("Handles.jl")

partial = (f::Function, y...)->(z...)->f(y..., z...)
∂(f::Function,y...) = partial(f, y...)
second(c) = c[2]

@HandleType HalfEdgeHandle
@HandleType VertexHandle
@HandleType FaceHandle
@HandleType EdgeHandle

using StaticArrays

const Vector3{T}=SVector{3,T}
const Vector3f=Vector3{Float32}
const Vector3d=Vector3{Float64}


""" a halfedge datastructure """
struct HalfEdge
  head     :: VertexHandle 
  next     :: HalfEdgeHandle
  opposite :: HalfEdgeHandle
  face     :: Union{FaceHandle,Nothing}
  HalfEdge( head_vert, next_edge, opp ) = new( head_vert, next_edge, opp, nothing )
  HalfEdge( head_vert, next_edge, opp, face ) = new( head_vert, next_edge, opp, face )
end

""" a manifold set of polygons """
struct Topology
  v2he    :: Vector{HalfEdgeHandle}
  face2he :: Vector{HalfEdgeHandle}
  he      :: Vector{HalfEdge}
  he2edge :: Vector{EdgeHandle}
end

const Mesh = Tuple{Topology, Vector{Vector3d}}
const Edge = Tuple{Int,Int}

HalfEdge(topo::Topology, heh::HalfEdgeHandle) = topo.he[heh]
""" dereference the HalfEdgeHandle to get the struct """
halfedge(topo::Topology, heh::HalfEdgeHandle)  = HalfEdge(topo, heh)
get(topo::Topology, heh::HalfEdgeHandle)  = HalfEdge(topo, heh)

""" get a halfedge for this face """
halfedge(topo::Topology, fh::FaceHandle) = topo.face2he[fh]

""" get a halfedge for this vertex """
halfedge(topo::Topology, vh::VertexHandle) = topo.v2he[vh]

"""
    halfedge( topo, vpair )

The oriented halfedge for the edge connecting two vertex indices/handles.
The halfedge points from first to second vertex in tuple vpair.
"""
halfedge(topo::Topology, e::T) where {T <: Tuple{VertexHandle, VertexHandle}} = 
  Iterators.filter(h->head(topo, h)==e[2], OneRing(topo, e[1])) |> first

halfedge(topo::Topology, e::T) where {V<:Integer, T<:Tuple{V,V}} = halfedge(topo, VertexHandle.(e))


HalfEdgeHandle(topo::Topology, vh::VertexHandle) = topo.v2he[vh]

""" return the VertexHandle for the vert pointed to by the given halfedge """
head(he::HalfEdge) = he.head
head(topo::Topology, heh::HalfEdgeHandle) = head(topo.he[heh])

""" return the VertexHandle for the vert pointed to by the given halfedge """
vertex(x::HalfEdge) = head(x)
vertex(topo::Topology, x) = head(topo, x)

vertices(topo::Topology) = topo.v2he
vertices(topo::Topology, c) = map(h->vertex(topo,h), c)
vertices(topo::Topology, heh::HalfEdgeHandle) = vertices(topo, Polygon(topo,heh))
vertices(topo::Topology, fh::FaceHandle) = vertices(topo, faces(topo)[fh])
verts(x...) = vertices(x...)

""" return the VertexHandle for the vert at the tail of the given halfedge """
tail(topo::Topology, he::HalfEdge) = head(topo,opposite(he))
tail(topo::Topology, heh::HalfEdgeHandle) = tail(topo,topo.he[heh])

""" return the next halfedge ccw on the same face """
function next(he::HalfEdge) 
  he.next
end

function next(topo::Topology, heh::HalfEdgeHandle) 
  next(topo.he[heh])
end

""" return halfedge pointing to this halfedge """
function prev(topo::Topology, heh::HalfEdgeHandle) 
  if isboundary(topo, heh)
    nothing
  else
    [ih for ih in Polygon(topo,heh)][end] 
  end
end

function opposite(he::HalfEdge) 
  he.opposite 
end

"""
    opposite(topo, halfedgehandle)

the Handle of the oppositely oriented halfedge sharing the same topological edge.   


        b
       /|\\
      /↑|↓\\ 
     / ↑|↓ \\
    ∘---a---∘

 ↑ and ↓ are opposite halfedges on the edge ab
"""
function opposite(topo::Topology, heh::HalfEdgeHandle) 
  opposite(topo.he[heh])
end

"iterate around a polygon"
struct Polygon
  topo::Topology
  starth::HalfEdgeHandle
end

Base.IteratorSize(::Type{Polygon}) = Base.SizeUnknown()
Base.eltype(::Type{Polygon}) = HalfEdgeHandle

""" iterator for a polygon """
Polygon(topo::Topology, fh::FaceHandle) = Polygon(topo, topo.face2he[fh])

Base.iterate( poly::Polygon ) = (poly.starth, next(poly.topo, poly.starth))

Base.iterate( poly::Polygon, heh::HalfEdgeHandle ) = 
  (heh == poly.starth) ?  nothing : (heh, next(poly.topo, heh))

function BoundaryLoop(topo::Topology, heh::HalfEdgeHandle) 
  if( isboundary(topo,heh) )
    Polygon(topo,heh)
  elseif( isboundary(topo,opposite(topo,heh)) )
    Polygon(topo,opposite(topo,heh))
  else
    nothing
  end
end

halfedges(topo) = HalfEdgeHandle.(1:length(topo.he))
halfedges(poly::Polygon) = [heh for heh in poly]

edge(topo::Topology, heh::HalfEdgeHandle) =  (tail(topo,heh), head(topo,heh))
"get an edge tuple from an EdgeHandle.  very slow"
edge(topo::Topology, eh::EdgeHandle) = edges(topo)[eh]
edgehandle(topo::Topology, heh::HalfEdgeHandle) = topo.he2edge[heh]


edges(topo::Topology) = Iterators.filter(ab->isless(ab...), Iterators.map(heh->edge(topo,heh),halfedges(topo))) |> collect
edges(topo::Topology, c) = map(h->edge(topo, h), c)
edges(topo::Topology, heh::HalfEdgeHandle) = edges(topo, Polygon(topo,heh))
edges(topo::Topology, fh::FaceHandle) = edges(topo, faces(topo)[fh])

mapedge(f, topo::Topology) = f.(Iterators.filter(ab->isless(ab...), Iterators.map(heh->edge(topo,heh),halfedges(topo))))


UniqueHalfEdges(topo::Topology) = Iterators.filter(heh->isless(edge(topo,heh)...), halfedges(topo))

"""
    UniqueHalfEdges(topo, flip_on_boundary)

Unique halfedges are halfedges whose ordering is in one-to-one correspondence with list of oriented edges.
we choose the edge (VertexHandle₁, VertexHandle₂) where VertexHandle₁ < VertexHandle₂.

i.e. map(heh->edge(topo,heh), UniqueHalfEdges(topo)) == edges(topo)

If flip_on_boundary is true then the one-to-one correspondence is not true when mapped to VertexHandle tuple.
In this case, the order of the UniqueHalfEdges still follows EdgeHandle order, but the orientation will be flipped at boundary.  This way you are guaranteed a valid face at each UniqueHalfEdge without having to check isboundary. 
"""
function UniqueHalfEdges(topo::Topology, flip_on_boundary)
  uhe = UniqueHalfEdges(topo)
  if flip_on_boundary
    Iterators.map(heh->isboundary(topo, heh) ? opposite(topo, heh) : heh, uhe)
  else
    uhe
  end
end

"""
unique halfedges not on boundary.
"""
unique_inside_halfedges(topo::Topology) = Iterators.filter(h->isboundary(topo, opposite(topo, h)) || (!isboundary(topo, h) && h  < opposite(topo, h)), halfedges(topo))

"unique indices for vertices"
vertex_indices(topo) = 1:length(topo.v2he)
"unique indices for faces"
face_indices(topo) = 1:length(topo.face2he)
"unique indices for edges"
edge_indices(topo) = 1:length(edges(topo)) # or length(topo.he)/2, but that is implemenation dependent

""" a boundary halfedge is connected to a single face """
isboundary(he::HalfEdge) =  he.face == nothing
isboundary(topo::Topology, heh::HalfEdgeHandle) =  isboundary(topo.he[heh])
isboundary(topo::Topology, heh::Nothing) =  true

""" a vertex is a boundary if it's starting halfedge is boundary """
# this assumes we always wind the halfedge to the boundary
# safer would be to iterate over the spokes and check if any halfedge is a boundary
isboundary(topo::T, vh::VertexHandle) where T = isboundary(topo, HalfEdgeHandle(topo, vh))

faces(topo::Topology) = topo.face2he

face(he::HalfEdge) =  he.face
face(topo::Topology, heh::HalfEdgeHandle) =  face(topo.he[heh])
face(topo::Topology, heh::Nothing) =  nothing

#!me can I reverse the order here so topology comes first?  Now that Handles are disticnt type I think I can
""" iterator for circulating clockwise around the one-ring neighbourhood of a vertex """
struct OneRing
  starth::HalfEdgeHandle
  topo::Topology
end

""" iterator to circulate around the one-ring neighbourhood of a vertex in cw direction """
OneRing(topo::Topology, vh::VertexHandle) = OneRing( HalfEdgeHandle(topo, vh), topo )

""" iterate around the one-ring of the vert in cw direction """
function Base.iterate(i::OneRing)
  heh = i.starth 
  (heh,next(i.topo,opposite(i.topo,heh)))
end

Base.iterate(topo::Topology, h::T) where T =  iterate(OneRing(topo,h))

# iterate around the one-ring of the vert in cw direction
function Base.iterate(i::OneRing, heh::HalfEdgeHandle)
  if heh == i.starth
    nothing
  else
    (heh,next(i.topo,opposite(i.topo,heh)))
  end
end

"""
  until( topo, he, v )

all halfedges from spoke a to b.  going ccw 
"""
function until( topo::Topology, he::HalfEdgeHandle, v::VertexHandle )
  r = Vector{HalfEdgeHandle}()
  push!( r, he )
  while( head(topo,he) != v )
    he = next( topo, he )
    push!( r, he )
  end
  r
end

""" 
  Rim(topo, i)

i is a VertexHandle

Iterates around rim of vertex ( outside edges of one-ring neighbourhood ). 
"""
function Rim(topo::Topology, i::VertexHandle) 
  spokeb = OneRing(topo,i)
  if isboundary(topo, first(spokeb))
    spokea = drop(spokeb, 1)
  else
    spokea = drop(cycle(spokeb), 1)
  end
  ( until(topo, next(topo, a), head(topo, b)) for (a,b) in zip(spokea, spokeb) ) |> flatten
end

#Base.length(i::OneRing) =  reduce((acc,_ )-> acc+1, i; init=0)
#!me compare performance of counting size first
Base.IteratorSize(::Type{OneRing}) = Base.SizeUnknown()

Base.eltype(i::OneRing) =  HalfEdgeHandle 

nhalfedges(topo) = length(topo.he)
nedges(topo) = nhalfedges(topo) >> 1
nvertices(topo) = length(vertices(topo))
nverts(topo) = nvertices(topo)
nfaces(topo) = length(faces(topo))


Topology( tris::Vector{Tuple{T,T,T}}, nVert::TT ) where {T<:Integer, TT<:Integer} = 
  Topology( map(t->(map(VertexHandle,t)...,),tris), nVert ) 

function Topology( poly::Vector{Vector{T}}, nVert::T; handle_bad_geo = true ) where {T<:Integer}
  nFace = length(poly)
  nHEdge = sum( map(length,poly) ) 
  topo_v2he, topo_he = (Vector{HalfEdgeHandle}(undef,nVert), Vector{HalfEdge}(undef,nHEdge)) 
  topo_face2he = Vector{HalfEdgeHandle}(undef,nFace)

  edgemap = Dict{Tuple{VertexHandle, VertexHandle}, HalfEdgeHandle}()
  heᵢ::HalfEdgeHandle = 1
  he2tail_vh = Vector{VertexHandle}(undef,nHEdge)

  for (i,fᵢ) in enumerate(poly) 

    nv = length(fᵢ)
    topo_face2he[i] = heᵢ
    # each half edge needs to link to next half edge. 
    # loop over edges and next halfedge handle, increment halfedge handle as we go
    heloop = zip(zip(drop(cycle(fᵢ), nv-1), fᵢ), take(drop(cycle(heᵢ:heᵢ+nv-1), 1), nv))
    for ((k,j),he_next) in heloop
      topo_he[heᵢ] = HalfEdge(j,he_next,he_next,FaceHandle(i)) 
      topo_v2he[j] = he_next
      edgemap[(k,j)] = heᵢ
      he2tail_vh[heᵢ] = k
      heᵢ += 1
    end

  end

  boundary_he = Vector{HalfEdge}(undef,0)
  topo_he2edge = Vector{EdgeHandle}(undef,0)
  prototopo = Topology(topo_v2he, topo_face2he, topo_he, topo_he2edge)
  topo_he = 
    map( zip(topo_he, 1:length(topo_he)) ) do (he,hehi) 
      headvert = head(he)
      tailvert = he2tail_vh[hehi]
      if haskey(edgemap, (headvert, tailvert)) 
        opp_heh = edgemap[(headvert, tailvert)] 
        HalfEdge(he.head, he.next, opp_heh, he.face ) 
      else 
        # create boundary halfedges. one on inside and one on outside.
        opposite_he = push!(boundary_he, 
                            HalfEdge(tailvert, HalfEdgeHandle(0), HalfEdgeHandle(hehi), nothing))
        opp_heh = HalfEdgeHandle(nHEdge + length(boundary_he))
        HalfEdge(he.head, he.next, opp_heh, he.face) 
      end
    end

  prototopo = Topology(topo_v2he, topo_face2he, vcat(topo_he, boundary_he), topo_he2edge)

  # iterate over boundaries.  set vertex map so boundary vertex tail is boundary he tail
  for i in 1:length(boundary_he)
    b_heh = HalfEdgeHandle(i+nHEdge)
    tailvert = head(prototopo, opposite(prototopo, b_heh))
    topo_v2he[tailvert] = b_heh
  end

  # now we can make a proper boundary halfedge
  boundary_he = 
    map( boundary_he ) do b_he
      # half edge from tail of the head vertex of this boundary he is next heh 
      # this *__>*
      #       ^  \  next
      #        \  v
      #          * 
      HalfEdge( head(b_he), topo_v2he[head(b_he)] ,opposite(b_he), nothing )
    end

  prototopo = Topology(topo_v2he, topo_face2he, vcat(topo_he, boundary_he), topo_he2edge)

  # make a map to EdgeHandles as this is quite useful for interfacing with linear algebra
  resize!(topo_he2edge, nhalfedges(prototopo))
  for (eh, hh) in enumerate(UniqueHalfEdges(prototopo))
    topo_he2edge[hh] = EdgeHandle(eh)
    topo_he2edge[opposite(prototopo, hh)] = EdgeHandle(eh)
  end

  Topology(topo_v2he, topo_face2he, vcat(topo_he, boundary_he), topo_he2edge)
end

function Topology( tris::Vector{T}, nVert::TT ) where {T<:Integer, TT<:Integer}
  trit = vec.(partition(tris, 3))
  if( length(trit[end]) < 3 )
    trit = trit[1:end-1]
  end
  Topology( trit, nVert ) 
end

"""
  Topology( P )

Represent a list of connected polygons P with halfedge datastructures.  P is either a flat list of triangles or a nested list of polygons with arbitrary vertex counts.  Polygons will be a list of vertex indices. 

For best results the polygons P should be manifold.

"""
Topology(polys) = Topology(polys,reduce(max,flatten(polys)))

# some useful geometry stuff
area(P,eab,ebc,eca) = norm((P[eab[2]] - P[eab[1]]) × (P[ebc[2]]-P[eab[1]]))*0.5

function area(P::T,eaz...) where {T<:AbstractArray}
  p₀ = P[eaz[1][2]]
  second = x->x[2]
  a = reduce(partition(second.(eaz[2:end]),2); init = 0.0) do acc,(i,j)
    acc+norm((P[i]-p₀) × (P[j]-p₀))
  end
  a*0.5
end

area(P::T,el) where {T<:AbstractArray} = area(P,el...)
area(topo::Topology, P, fh::FaceHandle) = area(P, edges(topo, fh)...) 
area(topo::Topology, P, heh::HalfEdgeHandle) = area(P, edges(topo, heh)...) 
area(topo::Topology, P, poly::Polygon) = area(topo, P, poly.starth)

function circumcenter_triangle(a::T,b::T,c::T) where T
  ac = c-a
  ab = b-a
  w = ab×ac

  u = (w×ab)*(ac⋅ac)
  v = (ac×w)*(ab⋅ab)
  x = (u+v)/(2.0 * w⋅w)

  x+a
end

function circumcenter_triangle(topo, P, h::HalfEdgeHandle)
  a = P[head(topo, h)]
  b = P[head(topo, next(topo, h))]

  if isboundary(topo, h) 
    return (a+b)*0.5
  end

  c = P[head(topo, next(topo, next(topo, h)))]
  circumcenter_triangle(a,b,c)
end

function dihedral_angle(topo::Topology, P, h)
  if isboundary(topo, h) || isboundary(topo, opposite(topo, h)) 
    return 0.0 
  end

  n1 = trinormal(topo, P, h)
  n2 = trinormal(topo, P, opposite(topo, h))
  w = (normalize∘-)(P[[edge(topo, h)...]]...)

  cosθ = n1⋅n2
  sinθ = (n1×n2)⋅w

  atan(sinθ, cosθ)
end 


"""
    inside_verts(topo, P, tris)

return the vertices on the inside of a set of overlapping tris
"""
function inside_verts(topo, P, overlaptris::Tuple{HalfEdgeHandle, HalfEdgeHandle})
#  W, height_cutoff = winding_number_cache(topo, P)

end

"""
    enclose(topo, P, boundary_faces)

return a list of all vertices enclosed by a contour of intersecting faces.  
If the mesh is not closed this may return all the vertices in the mesh. 
"""
function enclose(topo, P, boundary::Vector{Tuple{HalfEdgeHandle, HalfEdgeHandle}})
  #nPi = map(boundary_faces) do heh

  #end
end

#enclose(topo, P, bf::Vector{FaceHandle}) = enclose(topo, P, ∂(halfedge, topo).(bf))

"""
    floodfill(topo, P)

flood fill values at vertices where intersecting triangles create barrier.
"""
function floodfill(topo::Topology, P)
  # find intersections
  hits = collide_self(topo, P, triangle_edges);  
  # hit record example: (358, 6709, (true, Tuple{Int64,Tuple{Int64,Int64}}[(1, (2, 3)), (2, (1, 2))]))
  # (triangle1_index, triangle2_index, (true, [(which triangle, (edge vertices))...]))
  # means triangle 358 was pierced by edge with vertices polygons(topo)[358][[2,3]]
  # and   triangle 6709was pierced by edge with vertices polygons(topo)[6709][[1,2]]
  
  F = polygons(topo)
  n = normals(topo, P) 


  # this is not perfect, *--V--* hits 2 tris, can't determine which side it's on
  #!me need to use winding numbers
  for (tri1, tri2, (_, edges)) in hits
    for (whichtri, (a, b)) in edges
      if whichtri == 1
        facetri = tri1
        edgetri = tri2
      else
        facetri = tri2
        edgetri = tri1
      end
      Fedge = F[edgetri]
      Fface = F[facetri]
      side = n[facetri]⋅(P[Fedge[a]]-P[Fface[1]]) > 0.0 ? 1 : -1
      hitside[Fedge[a]] = side
      side = n[facetri]⋅(P[Fedge[b]]-P[Fface[1]]) > 0.0 ? 1 : -1
      hitside[Fedge[b]] = side
    end
  end

  # isolate islands

end

"""
    polygons(topo)

extract polygons as arrays of vertex indices 
"""
polygons(topo) = map(f->vertices(topo,Polygon(topo,f)), (faces(topo)))

minussigned( ::Val{false} ) = 1.0
minussigned( ::Val{true} ) = -1.0

"""
    incidence(topo, AHandleType, BHandleType, oriented = false)

build incidence matrix mapping between two elements
|b|x|a| matrix. i.e. incidence of {vertex} and {face} has size(vertex) columns and size(face) rows

If oriented then matrix element will be negative if elements are oriented in opposite sense.
For a vertex-edge this means positive if the vertex is at the head of the edge.
For edge-face this means positive if direction of circulation of the edge direction is ccw relative to face
"""
function incidence(topo, a::Type{VertexHandle}, b::Type{EdgeHandle}, oriented = false)
  R = Int[]
  C = Int[]
  for (eᵢⱼ,(vᵢ,vⱼ)) in enumerate(edges(topo))
    push!(C, vᵢ)
    push!(C, vⱼ)
    push!(R, eᵢⱼ)
    push!(R, eᵢⱼ)
  end
 
  if oriented
    V = map(i->(i%2 ==1) ? -1.0 : 1.0, 1:length(C))
  else
    V = ones(Float64, length(C))
  end

  sparse(R,C,V)
end

function incidence(topo, a::Type{EdgeHandle}, b::Type{FaceHandle}, oriented = false)
  R = Int[]
  C = Int[]
  V = Float64[]
  for (eᵢⱼ,heᵢⱼ) in enumerate(UniqueHalfEdges(topo))
    f = face(topo, heᵢⱼ)
    if f != nothing
      push!(R, f)
      push!(C, eᵢⱼ)
      push!(V, 1.0)
    end
    f = face(topo, opposite(topo,heᵢⱼ))
    if f != nothing
      push!(R, f)
      push!(C, eᵢⱼ)
      push!(V, minussigned(Val(oriented)))
    end
  end
  
  sparse(R,C,V)
end

hashedge(topo::Topology, a, b) = Int((b << (leading_zeros(0)-leading_zeros(nedges(topo)) + 1)) + a)

hashedge(topo::Topology, (a,b)) = hashedge(topo, a, b)

function unhashedge(topo::Topology, ab::Int)
  bshift = (leading_zeros(0)-leading_zeros(nedges(topo)) + 1)
  b = ab >> bshift
  a = ab - (b << bshift)
  (a,b)
end

Base.sort((a,b)::T) where {V,T<:Tuple{V,V}} = a>b ? (b,a) : (a,b)

"""
  IncidentEdges( topo, HandleType )

Get list of incident edges as EdgeHandle[]
Edges will be ordered correctly with winding order of face.
"""
function IncidentEdges(topo::Topology, ::Type{FaceHandle})
  #!me can rewrite this to using he2edge map

  # first place hash of edges in correctly ordered positions
  incidente = Vector{Vector{EdgeHandle}}() 
  sizehint!(incidente, nfaces(topo))
  for shf in faces(topo)
    ev = Vector{EdgeHandle}()
    sizehint!(ev,3)
    nhf = shf
    while true
      ab = edge(topo, nhf)
      push!(ev, EdgeHandle(hashedge(topo, sort(ab))))
      nhf = next(topo, nhf)
      nhf == shf && break
    end
    push!(incidente, ev)
  end

  # find and replace hashed edges with unique edge handles
  for (ie, heh) in enumerate(UniqueHalfEdges(topo))
    ehsh = hashedge(topo, sort(edge(topo, heh)))
    if !isboundary(topo, heh)
      replace!(incidente[face(topo, heh)], ehsh=>ie; count=1)
    end
    oheh = opposite(topo, heh)
    if !isboundary(topo, oheh)
      replace!(incidente[face(topo, oheh)], ehsh=>ie; count=1)
    end
  end

  incidente
end

"""
  IncidentHalfEdges( topo, HandleType )

Get list of incident halfedges as HalfEdgeHandle[]
Halfedges will be ordered correctly with winding order of face.
"""
function IncidentHalfEdges(topo::Topology, ::Type{FaceHandle})
  # first place hash of edges in correctly ordered positions
  incidente = Vector{Vector{HalfEdgeHandle}}() 
  sizehint!(incidente, nfaces(topo))
  for shf in faces(topo)
    ev = Vector{HalfEdgeHandle}()
    sizehint!(ev,3)
    nhf = shf
    while true
      ab = edge(topo, nhf)
      push!(ev, HalfEdgeHandle(hashedge(topo, sort(ab))))
      nhf = next(topo, nhf)
      nhf == shf && break
    end
    push!(incidente, ev)
  end

  # find and replace hashed edges with unique halfedge handles
  for heh in UniqueHalfEdges(topo)
    ehsh = hashedge(topo, sort(edge(topo, heh)))
    if !isboundary(topo, heh)
      replace!(incidente[face(topo, heh)], ehsh=>heh; count=1)
    end
    oheh = opposite(topo, heh)
    if !isboundary(topo, oheh)
      replace!(incidente[face(topo, oheh)], ehsh=>heh; count=1)
    end
  end

  incidente
end

"""
GeometryTypes mesh to halfedge and extracted tris, verts
"""
function geometry2hemesh( cat )
  tris = map(t->[map(i->Int(i)+1,t)...],cat.faces)
  P = Vector3d.(cat.vertices)
  hemesh = Topology(tris)
  (hemesh,tris,P)
end

"""
  facelist(topo)

return faces as arrays of integers
"""
facelist(topo::Topology) = (x->Int.(vertices(topo,Polygon(topo,x)))).(faces(topo))


"""
  trinormal(topo, P, h)

normal for triangle associated with HalfEdgeHandle h 
"""
function trinormal(topo::Topology, P, h::HalfEdgeHandle)
  h = isboundary(topo, h) ? opposite(topo, h) : h
  pa = P[head(topo, h)]; h = next(topo, h)
  pb = P[head(topo, h)]; h = next(topo, h)
  pc = P[head(topo, h)]
  n = (pb-pa)×(pc-pa)
  nn = n⋅n
  if isapprox(nn, 0, atol=eps())
    normalize(ones(pa))
  else
    n/sqrt(nn)
  end
end

"""
    normals(topo, P)

triangle normals for all faces
"""
normals(topo::Topology, P) = map( fheh->trinormal(topo, P, fheh), faces(topo)) 


"""
    weightednormal(topo, P, heh)

area weighted normal of a triangle
"""
function weightednormal(mesh, P, heh::HalfEdgeHandle)
  p₀ = P[head(mesh, heh)]
  heh₁ = next(mesh, heh)
  p₁ = P[head(mesh, heh₁)]
  p₂ = P[head(mesh, next(mesh, heh₁))]

  0.5*(p₁-p₀)×(p₂-p₀)
end

"""
    vertexnormal(topo, p, ring)

normal of a vertex given it's OneRing
"""
function vertexnormal( mesh, P, ring::OneRing, normal = weightednormal )
  iring = Iterators.filter(heh->!isboundary(mesh,heh),ring)
  normalize(mapreduce(heh->normal(mesh,P,heh),+,iring))
end

"""
    vertexnormals(topo, P, weighted=true)

vertex normals for all vertices in mesh, optionally area weighted.    
"""
function vertexnormals( mesh, P, weighted = true )
  rings = (heh->OneRing(heh,mesh)).(vertices(mesh))
  if weighted
    n = ring->vertexnormal( mesh, P, ring, weightednormal )
  else
    n = ring->vertexnormal( mesh, P, ring, normalize∘weightednormal )
  end
  n.(rings)
end

"""
  angles(topo, P, poly)

Interior angles for polygon as ccw rotation around tail of each halfedge
"""
function angles(topo::Topology, P, poly::Polygon) 
  n = length(poly|>collect)
  hv = (head(topo,h) for h in cycle(poly))
  abc = zip(drop(hv,n-1), hv, drop(hv,n-2))
  map(take(abc,n)) do (a,b,c)
    acos(normalize(P[b]-P[a])⋅normalize(P[c]-P[a]))
  end
end

"""
  angle(topo, P, heh)

Interior angle across from the halfedge
"""
function angle(topo, P, h::HalfEdgeHandle)
  poly = Polygon(topo, h)
  n = length(poly|>collect)
  nop = n >> 1
  hv = (head(topo,h) for h in cycle(poly))
  a,b,c = P[[first(drop(hv,nop)), first(drop(hv,nop-1)), first(drop(hv,nop+1))]]
  acos(normalize(b-a)⋅normalize(c-a))
end

"""
    orientation( (a,b,c), (u,v) )

returns 1.0 if the edge u,v has same winding as a,b,c otherwise -1.0
"""
function orientation( a::H ,b::H, c::H, u::H, v::H ) where {H<:Union{VertexHandle, HalfEdgeHandle}}
  (u==a && v==b) || ((u==b) && (v==c)) || ((u==c) && (v==a)) ? 1.0 : -1.0
end

orientation( abc::FT, uv::ET ) where {FT<:Union{Vector, Tuple}, ET<:Union{Vector, Tuple}} = 
  orientation(abc..., uv...)

"""
  improve_mesh(topo, P; smoothness)

Improve the quality of the mesh.
Particularly so discrete differential geometry operators yield well conditioned matrices 
Note: input topology will not longer match.  i.e. FaceHandle(i) will not necessarily return same vert tuple

Currently element counts remain the same.  The only operation performed is to flip edges between two triangles so that all angles are closer to π/2


Returns the number of edges flipped
"""
#==
#!me  not really working
If smoothness is set to something other than 0, edge flips will not occur if the volume of the tetrahedron formed from the 4 vertices involved as a ratio against triangle areas is less than the smoothness.
In other words, smoothness helps preserve the silhouette of the mesh. 

==#
function improve_mesh(topo::Topology, P) #; smoothness=0.0)
  # reverse map into Topology::v2he and Topology::face2he data
  he2v = Vector{Union{Missing,Int}}(missing, length(topo.he))
  he2face = Vector{Union{Missing,Int}}(missing,length(topo.he))

  for (i,heᵢ) in enumerate(topo.v2he)
    he2v[heᵢ] = i 
  end

  for (i,heᵢ) in enumerate(topo.face2he)
    he2face[heᵢ] = i 
  end

  for heᵢ in halfedges(topo)
    if ismissing(he2face[heᵢ])
      for hpᵢ in Polygon(topo, heᵢ)
        if !ismissing(he2face[hpᵢ])
          he2face[heᵢ] = he2face[hpᵢ]
          break;
        end
      end
    end
  end

  hasflipped = falses(length(topo.face2he))

  for heᵢ in unique_inside_halfedges(topo)
    bscale = 0.0
    α = 0.0
    β = 0.0
    skipα = false
    skipβ = false
    oheᵢ = opposite(topo, heᵢ)

    if hasflipped[he2face[heᵢ]] || 
      (!isboundary(topo, oheᵢ) && hasflipped[he2face[oheᵢ]])
      continue
    end

    if isboundary(topo, oheᵢ) || length(Polygon(topo, he2face[oheᵢ])|>collect) != 3
      skipβ = true
    else
      β = angle(topo, P, oheᵢ)
      polyβ = he2face[oheᵢ]
      heflip = oheᵢ
      bscale = 0.5
    end
    if isboundary(topo, heᵢ) || length(Polygon(topo, he2face[heᵢ])|>collect) != 3
      skipα = true
    else
      α = angle(topo, P, heᵢ)
      polyα = he2face[heᵢ]
      heflip = heᵢ
      bscale += 0.5
    end


    if α + β > π*bscale

      if !skipα && !skipβ
        #==
            N
           /α\
          / ↳ \
         W=⃭=⃭=⃭=⃭=⃭E
          \ ↰ /
           \β/
            S

            to

            N
           /|\
          / | \
         W ↰| ↰E 
          \α|β/
           \|/
            S
        ==#
        hmidα = heᵢ; midα = topo.he[heᵢ] 
        hmidβ = oheᵢ; midβ = topo.he[oheᵢ] 
        hNE = midα.next; NE = topo.he[hNE]
        hSW = midβ.next; SW = topo.he[hSW]
        hNW = NE.next; NW = topo.he[hNW]
        hSE = SW.next; SE = topo.he[hSE]

        N = NE.head; W = NW.head; S = SW.head; E = SE.head

        # skip if the silhouette will be made significantly less smooth
        #==
        if smoothness > 0.0
          # check ratio of volume of tetrahedron vs avg area
          n = P[N]; s = P[S]; e = P[E]; w = P[W]
          se = e-s; sw = w-s; 
          base = se×sw
          h = normalize(base)⋅(n-s)
          avgl = (norm(n-s) + norm(w-e))*0.5
          
          if h/avgl > 1.0/smoothness
            continue
          end
        end
        ==#

        # new midα. overwrite midα
        topo.he[heᵢ] = HalfEdge(N, hNW, hmidβ, midα.face)
        # new midβ. overwrite midβ
        topo.he[oheᵢ] = HalfEdge(S, hSE, hmidα, midβ.face)

        # new NE
        topo.he[hNE] = HalfEdge(N, hmidβ, NE.opposite, midβ.face)
        # new NW
        topo.he[hNW] = HalfEdge(W, hSW, NW.opposite, midα.face) 
        # new SW
        topo.he[hSW] = HalfEdge(S, hmidα, SW.opposite, midα.face) 
        # new SE
        topo.he[hSE] = HalfEdge(E, hNE, SE.opposite, midβ.face) 

        # fix facelist
        topo.face2he[he2face[hmidα]] = hmidα
        topo.face2he[he2face[hmidβ]] = hmidβ

        # fix vertexlist.  keep in mind iteration from boundary, so we can't pick arbitrary he
        if topo.v2he[E] == hmidβ
          topo.v2he[E] = hNE
        end
        if topo.v2he[W] == hmidα
          topo.v2he[W] = hSW
        end

        hasflipped[polyα] = true
        hasflipped[polyβ] = true

      else
        # ignore boundary for now as it requires adding a new point and face
      end

    end
  end
  sum(hasflipped)
end

#==
struct Manifold
  soup2manifold::Vector{VertexHandle}
  manifold2soup::Vector{VertexHandle}
  manifoldpoly::Vector{Vector{VertexHandle}}

  orphan::Vector{VertexHandle}
end

function make_manifold( P, poly )
  nVert = length(P)

  # weld overlapping verts
  # AABBTree()

  # fix any fins or other non-manifold topology 
  #

  # isolate islands and keep the biggest one?
    
  orphanv = VertexHandle[]
  floaters = setdiff(1:nVert,reduce(vcat,poly))
  for i in floaters
    push!(orphanv, i)
  end


end
==#


"""
  boundary_verts(topo)

find loops of boundary verts
"""
function boundary_verts(topo) 
  bh = filter(∂(isboundary, topo), halfedges(topo))

  bvl = Vector{VertexHandle}[]

  while !isempty(bh)
    bl = BoundaryLoop(topo, bh |> first)
    bh = setdiff(bh, bl)
    push!(bvl, map(∂(tail, topo), bl))
  end

  bvl
end

boundary_vertices(x) = boundary_verts(x)

"""
  boundary_interior(topo, boundary)

the set of vertices connected to the boundary verts but not part of the boundary
"""
function boundary_interior(mesh, boundary)
  allverts = mapreduce(x->OneRing(mesh, x) |> ∂(∂(mapreduce,(x->[x...])∘∂(edge,mesh)),vcat), vcat , boundary)
  setdiff(allverts,boundary)
end

boundary_interior(topo) = map(b->boundary_interior(topo, b), boundary_verts(topo))

# single airty versions

for fn in (:next, :head, :tail, :prev, :opposite, :isboundary, :edge, :halfedge)
  @eval $fn(h) = topo->$fn(topo, h)
  @eval $fn(f::F) where {F<:Function} = (topo,h)->$fn(topo,f(topo,h))
end

"""
parse an .obj file format mesh into points and faces.
"""
function parseobj(instream)
  # of form "f v1 v2 v3 ....""
  process_face(words::Vector{S}) where {S <: AbstractString} = (words,) # make same format as the others
  # of form "f v1//vn1 v2//vn2 v3//vn3 ..."
  process_face_normal(words::Vector{S}) where {S <: AbstractString} = mapreduce(first,vcat,split.(words, "//"))
  # of form "f v1/vt1 v2/vt2 v3/vt3 ..." or of form "f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3 ...."
  process_face_uv_or_normal(words::Vector{S}) where {S <: AbstractString} = mapreduce(first,vcat,split.(words, Ref('/')))
  process_face_uv(words::Vector{S}) where {S <: AbstractString} = mapreduce(second,vcat,split.(words, Ref('/')))

  # this is cribbed from https://github.com/JuliaIO/MeshIO.jl/blob/master/src/io/obj.jl
  #  normals, uvs etc can be parsed using that package
  v = Vector3d[]
  uv = Vector{Float64}[]
  f = Vector{Int}[]
  f_uv = Vector{Int}[]
  for line in filter(ln->!(startswith(ln,"#") || 
                           isempty(ln) || 
                           all(iscntrl,ln)), 
                     (strip∘chomp).(eachline(instream)))
    # assuming one element per line
    words = split(line)
    command = popfirst!(words)

    if command == "v"
      push!(v, Vector3d(parse.(Float64, words)))
    elseif command == "vt"
      push!(uv, parse.(Float64, words))
    elseif command == "f"
      if any(x->occursin("//", x), words)
        fs = process_face_normal(words)
        push!(f,parse.(Int,fs))
      elseif any(x->occursin("/", x), words)
        fs = process_face_uv_or_normal(words)
        push!(f,parse.(Int,fs))
        fs_uv = process_face_uv(words)
        push!(f_uv,parse.(Int,fs_uv))
      else
        fs = words
        push!(f,parse.(Int,fs))
      end
    end

  end

  # obj uses base 1 indexing
  mesh = Dict("P"=>v, "Poly"=>f)
  if length(f_uv) > 0
    mesh["PolyUV"] = f_uv
    mesh["UV"] = uv
  end

  mesh
end

"""
    loadmesh("pathto.obj")

Load a Mesh from an obj file.  Topology and points.    
"""
function loadmesh(filename::AbstractString)
  if split( filename, "." )[end] == "obj"
    local topo, P
    open(filename, "r") do f
      rp = parseobj(f)
      topo = Topology(rp["Poly"]); 
      P = rp["P"];
    end
    (topo, P) 
  else
    :FileTypeNotSupported 
  end
end

#==========  Collision Detection ============#

include("BoundingVolumeTrees.jl")
include("Collision.jl")
include("IslandFind.jl")

using .BoundingVolumeTrees
using .Collision
using .IslandFind


export 
BVH,
Collider,
collide_self,
query_aabb,
update_collider

const BVH = BoundingVolumeTrees

struct Collider{T,F}
  bvh::BVH.DAABBNode{T,F}
  mesh::BVH.Triangles{T,F}
end

"""
    Collider(topo, P)

precompute collision detection accelerators
"""
function Collider(topo::Topology, P::V) where {T, PT<:AbstractVector{T}, V<:Vector{PT}}
  mesh = BVH.Triangles((collect∘flatten)(facelist(topo)),P)
  bvh = BVH.createBVH(mesh)
  Collider{Int64,T}(bvh, mesh)
end

"""
    collide_self(topo, P, collidef)

return a list of faces which are intersecting in the mesh.
provide the face vs face collision method in collidef

collidef should be a function that accepts 6 points as arguments, which are the points of the two triangle.
It should return a tuple where the first element is a boolean value indicating if there was a collision or not
"""
collide_self(topo::Topology, P, collidef::F = Collision.triangle_triangle) where {F<:Function} = collide_self(Collider(topo, P), collidef) 

"""
    collide_self(collider, collidef)

return a list of faces which are intersecting in the mesh.
provide the face vs face collision method in collidef
"""
function collide_self(collider::Collider, collidef::F = Collision.triangle_triangle) where{F<:Function}
  BVH.selfintersect(collider.bvh, collider.mesh, collidef)
end

"""
    query_aabb(collider, aabb, hits)

query the collision detection structure for a list of faces overlapping the given axis aligned bounding box
new hits are pushed onto the passed in hits vector
"""
query_aabb(c::Collider{T,F}, aabb::BVH.AABB{F}, hits::Vector{T}) where {T,F} = BVH.query(c.bvh, aabb, hits) 

function query_aabb(c::Collider{T,F}, aabb::BVH.AABB{F}) where {T,F}
  hits = Vector{T}()
  BVH.query(c.bvh, aabb, hits)
  hits
end

"""
    update_collider(collider)

update the collision detection structures in collider to reflect any change in mesh positions
"""
update_collider( collider::C ) where C<:Collider = BVH.updateBVH(collider.bvh, collider.mesh)

include("WindingNumbers.jl")

end # module
