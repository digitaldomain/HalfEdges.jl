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

using Test
using HalfEdges
using HalfEdges: halfedges, vertices, vertex, edges, edge, halfedge, face, head, tail, opposite, next
using HalfEdges: orientation, polygons
using Serialization
using StaticArrays

∂(f::Function, y...) = (z...)->f(y...,z...)
second(c) = c[2]
he_ = HalfEdges

@testset "one tri" begin
  topo = HalfEdges.Topology([1,2,3])
  @test 1 == 1  # just make sure it didn't fail out
  @test vertex(HalfEdge(topo, opposite(topo, HalfEdgeHandle(topo, VertexHandle(1))))) == 1
  @test tail(topo, vertices(topo)[1])  == VertexHandle(1)
  @test he_.mapedge( i->i, topo)[1] == edges(topo)[1]
  @test he_.vertex_indices(topo) == 1:nverts(topo)
  @test he_.face_indices(topo) == 1:nfaces(topo)
  @test sort(vertices(topo, HalfEdgeHandle(1))) == sort(vertices(topo, FaceHandle(1)))
  @test length(halfedges(Polygon(topo, 1))) == 3
  @test eltype(OneRing(topo,VertexHandle(1))) == HalfEdgeHandle
  @test HalfEdges.verts(topo) == vertices(topo)
  @test HalfEdges.edge_indices(topo) == 1:nedges(topo)

  @test isboundary(topo, nothing)
  @test face(topo, nothing) == nothing

  P = [[0.0,0.0,0.0], [1.0,0.0,0.0], [0.0,1.0,0.0]]
  hface = opposite(topo, HalfEdgeHandle(topo, VertexHandle(1)))
  @test isboundary(topo, hface) == false
  @test he_.circumcenter_triangle(topo, P, hface) == [0.5,0.5,0.0]
  @test he_.normals(topo, P) == [[0.0,0.0,1.0]]
  @test facelist(topo) == [[1,2,3]]
  @test edge(topo, EdgeHandle(1)) == (1,2)
  @test sort(he_.angles(topo, P, Polygon(topo, 1))) ≈ [π/4, π/4, π/2]

  # HalfEdgeHandle for vertex 3 is the most counterclockwise handle with 3 as it's tail
  @test isboundary(topo, vertices(topo)[3])
  @test edge(topo, vertices(topo)[3]) == (VertexHandle(3), VertexHandle(2))
  @test he_.angle(topo, P, vertices(topo)[3]) ≈ π/2
  @test he_.edgehandle(topo, first(UniqueHalfEdges(topo))) ==
            he_.edgehandle(topo, opposite(topo, first(UniqueHalfEdges(topo))))

  @test he_.area(topo, P, FaceHandle(1)) == 0.5
  @test isboundary(topo, VertexHandle(1)) 
end


@testset "bounded ring" begin
  topo = Topology([[1,2,3],[3,2,4],[4,5,3]],5)

  for i in 1:length(topo.v2he)
    heh = topo.v2he[i]
    @test heh == HalfEdges.opposite(topo, HalfEdges.opposite(topo, heh))
  end

  @test length(collect(HalfEdges.OneRing(topo,VertexHandle(1)))) == 2
  @test sort(OneRingVerts(topo,VertexHandle(1))) == [2, 3]
  @test length(collect(HalfEdges.OneRing(topo,VertexHandle(3)))) == 4
  @test length(collect(HalfEdges.OneRing(topo,VertexHandle(4)))) == 3
  ring3 = [HalfEdges.edge(topo,ih) for ih in HalfEdges.OneRing(topo,VertexHandle(3))]
  @test reduce( (acc,iv)->acc && iv == VertexHandle(3), [HalfEdges.tail(topo,ih) for ih in HalfEdges.OneRing(topo,VertexHandle(3))][1:end], init=true)
  # testing assumption for BoundaryLoop test, not a "real" test ( just assuming hedge 1 is boundary)
  @test HalfEdges.isboundary(topo, HalfEdgeHandle(1)) || 
        HalfEdges.isboundary(topo, HalfEdges.opposite(topo,HalfEdgeHandle(1)))
  @test HalfEdges.BoundaryLoop(topo,HalfEdgeHandle(1)) |> collect |> length == 5
  @test nhalfedges(topo) == length(edges(topo))*2
  alledge = edges(topo); allface = polygons(topo);
  @test orientation(allface[1], alledge[1]) == 1.0
  @test orientation(allface[1], (VertexHandle(2), VertexHandle(3))) == 1.0
  @test orientation(allface[2], (VertexHandle(2), VertexHandle(3))) == -1.0
  @test edges(topo)[IncidentEdges(topo, FaceHandle)[1]] == map(sort, edges(topo, FaceHandle(1)))
  @test edges(topo)[IncidentEdges(topo, FaceHandle)[2]] == map(sort, edges(topo, FaceHandle(2)))
  @test edges(topo)[IncidentEdges(topo, FaceHandle)[3]] == map(sort, edges(topo, FaceHandle(3)))
  @test map(heh->edge(topo,heh), IncidentHalfEdges(topo, FaceHandle)[3]) == map(sort, edges(topo, FaceHandle(3)))
end


@testset "unbroken ring" begin
  topo = Topology([[1,2,3],[1,3,4],[1,4,2]],4)
  P = [[0.0,0.0,0.0],[1.0,-0.25,0.0],[0.0,1.0,0.0],[-1.0,-0.25,0.0]]

  @test length(collect(HalfEdges.OneRing(topo,VertexHandle(1)))) == 3
  @test length(collect(HalfEdges.OneRing(topo,VertexHandle(4)))) == 3
  @test reduce( (acc,iv)->acc && iv == VertexHandle(1), [HalfEdges.tail(topo,ih) for ih in HalfEdges.OneRing(topo,VertexHandle(1))][1:end], init=true)
  @test nhalfedges(topo) == length(edges(topo))*2

  FE = he_.incidence(topo, EdgeHandle, FaceHandle,true)
  EV = he_.incidence(topo, VertexHandle, EdgeHandle,true)
  @test iszero(FE*EV) 
  @test length(vertexnormals(topo, P)) == nvertices(topo)
  @test vertexnormals(topo,P)[1] == [0.0, 0.0, 1.0]
  @test length(normals(topo, P)) == nfaces(topo) 
  @test normals(topo,P)[1] == [0.0, 0.0, 1.0]
end

@testset "3x3grid" begin
  cell(i,j) = map(v->v+(j-1)*4, [i,i+4,i+1,i+1,i+4,i+5])
  topo = Topology(Iterators.flatten([cell(1,1), cell(2,1), cell(3,1), 
                                     cell(1,2), cell(2,2), cell(3,2), 
                                     cell(1,3), cell(2,3), cell(3,3)])|>collect)
  @test sort(first(he_.boundary_interior(topo))) == [6,7,10,11]

end

@testset "box48" begin
  f = open(string(Base.@__DIR__, "/box48.spork"), "r")
  mesh = deserialize(f)
  close(f)

  hemesh = Topology(mesh["Poly"], length(mesh["P"]))
  P = mesh["P"]

  @test sum(map(i->he_.area(hemesh,P,i), he_.faces(hemesh))) == 6.0
  @test he_.nedges(hemesh) == he_.edges(hemesh) |> length
  @test he_.nhalfedges(hemesh)/2 == he_.nedges(hemesh) == he_.length(he_.edges(hemesh))
  @test he_.edges(hemesh) |> collect == ∂(edge,hemesh).(he_.UniqueHalfEdges(hemesh)) 
end

@testset "incidence" begin
  topo = Topology([[1,2,3],[2,4,3]])

  @test edges(topo)[IncidentEdges(topo, FaceHandle)[1]] == map(sort, edges(topo, FaceHandle(1)))
  @test edges(topo)[IncidentEdges(topo, FaceHandle)[2]] == map(sort, edges(topo, FaceHandle(2)))

  P = [[0.0,0.0,0.0], [1.0,0.0,0.0], [0.0,1.0,0.0], [0.5, 0.5,1.0]]
  @test sort(map(h->he_.dihedral_angle(topo, P, h), UniqueHalfEdges(topo)))[[1,end]] ≈ [0.0,π/2]
  # not really working so well, but we will test it anyways
  P = [[0.25,0.25,0.0], [1.0,0.0,0.0], [0.0,1.0,0.0], [0.75, 0.75,0.0]]
  @test he_.improve_mesh(topo, P) > 0
end

@testset "readme" begin
  topo = Topology([[1,2,3],[1,4,2]])

  @test edges(topo) |> collect == ∂(edge, topo).(UniqueHalfEdges(topo))

  # find the only halfedge pointing from vertex 1 to 2
  h = Iterators.filter(h->head(topo, h)==VertexHandle(2), OneRing(topo, VertexHandle(1))) |> first
  @test h == halfedge(topo,(1,2))
  @test tail(topo, h) == VertexHandle(1)
  @test head(topo, h) == VertexHandle(2)
  @test vertex(topo, h) == VertexHandle(2)

  h = opposite(topo, h)
  @test h == halfedge(topo, (2,1))

  topo = Topology([[1,2,3],[1,4,2],[1,3,5,4],[2,4,6]])
  @test length(edges(topo, OneRing(topo, VertexHandle(1)))) == 3
  @test length(edges(topo, Rim(topo, VertexHandle(1)))) == 4
  @test length(edges(topo, BoundaryLoop(topo, first(OneRing(topo, VertexHandle(2)))))) == 5
  @test length(edges(topo, Polygon(topo, FaceHandle(1)))) == 3
  @test length(edges(topo, Polygon(topo, FaceHandle(3)))) == 4


  v6 = VertexHandle(6)
  @test head(topo, next(topo, opposite(topo, next(topo, next(topo, opposite(topo, next(topo, HalfEdgeHandle(1)))))))) == v6
  pivot = opposite(next(next))

  @test head(next)(topo, pivot(topo, opposite(next)(topo, HalfEdgeHandle(1)))) == v6

  @test (halfedge |> next |> opposite |> next |> next |> opposite |> next |> head)(topo, (3,1)) == v6 
end

@testset "Degenerate" begin
  topo = Topology([2,3,4])
  P = map(i->rand(SVector{3,Float64}), 1:4)
  @test first(vertexnormals(topo, P)) == zero(SVector{3,Float64})
  @test isempty(floodfill(topo, P)) == false
end

@HandleType Foo
@HandleType Bar
@HandleType Baz

hey(x::Integer) = x+1
hey(x::Foo) = string("foo", x)
hey(x::Bar) = string("bar", x)

@testset "Handles" begin
  @test Foo(1) == 1

  a = sort(map( f->f+1 + Foo(1), (2 .* Foo.([10:-1:-2...]))))
  @test eltype(a) == Foo
  @test a[1] == -2 == Foo(-2)
  @test a[2] < a[3]
  @test a[end] == 22

  @test hey(Foo(2)) == "foo2"
  @test hey(Bar(3)) == "bar3"
  @test hey(Baz(4)) == 5

end

function cube()
  P = [0 0 0
       1 0 0
       0 0 1
       1 0 1
       1 1 1
       0 1 1
       1 1 0
       0 1 0
      ]
  P = SVector{3}.(eachrow(Float64.(P)))
  (Topology([1,2,3, 2,4,3, 3,4,6, 4,5,6, 2,5,4, 2,7,5, 7,8,5, 5,8,6, 1,3,6, 1,6,8, 1,8,2, 2,8,7]), P)
end

@testset "collide" begin
  topo,P = cube()
  col = Collider(topo, P)
  tophalf = query_aabb(col, HalfEdges.BVH.AABB(SVector{3}(0.,0.,0.5), SVector{3}(1.0,1.0,1.0)))  
  pre = []
  inorder = []
  post = []
  HalfEdges.BVH.visitdf(col.bvh, x->push!(pre, x), x->push!(post, x), x->push!(inorder, x))
  @test length(pre) == length(post) == length(inorder) > 0

  # should contain all triangles other than bottom two
  @test length(tophalf) == nfaces(topo)-2
  @test isempty(tophalf ∩ [11, 12])   

  #pull point into self intersection
  col.mesh.P[1] = SVector{3}(0.5,0.5,2.0)
  update_collider(col)
  #foo(a1,a2,a3, b1,b2,b3) = a2[1] > 1.0
  #@test collide_self(col, foo) == collide_self(topo, col.mesh.P, foo)
  @test collide_self(col) == collide_self(topo, col.mesh.P)

  #==
      2---4
     /|\  |\
    5 | \ | 6    then move 5 and 6 so they intersect
     \|  \|/
      1---3
  ==#
  F = [[1,3,2],[2,3,4],[2,5,1],[3,6,4]]
  topo = Topology(F)
  P = [0 0 0
       0 1 0
       1 0 0
       1 1 0
       1 0.5 1
       0 0.5 1
      ]
  P = SVector{3}.(eachrow(Float64.(P)))
  #foo(a1,a2,a3, b1,b2,b3) = a1[2] > 0.0   # our two triangle with z > 0 are in collision
  #hits = collide_self(Collider(topo, P), foo)
  hits = collide_self(Collider(topo, P))
  @test length(hits) == 1

  # want to get lots of hits to test when hits array is resized
  topo = Topology(repeat(F, 100))
  #hits = collide_self(Collider(topo, P), foo)
  hits = collide_self(Collider(topo, P))
  @test length(hits) > 1
end

@testset "winding numbers" begin
  topo, P = cube()

  Fi = FaceHandle.(1:nfaces(topo))
  centre = sum(P)/length(P)
  @test sum(map(fh->he_.solid_angle(topo, P, he_.halfedge(topo, fh), centre), Fi)) ≈ 4.0*π
  @test sum(map(fh->he_.solid_angle(topo, P, he_.halfedge(topo, fh), P[1]), Fi)) ≈ π/2.0

  # outside, accuracy seems lower for some reason
  @test abs(sum(map(fh->he_.solid_angle(topo, P, he_.halfedge(topo, fh), P[1]-centre), Fi))) < 1e6 

  cache = winding_number_cache(topo, P, 1)
  #cache3 = winding_number_cache(topo, P, 1, 3)

  @test winding_number(topo, P, SVector{3}(-1.0,-1.0,-1.0), cache) |> round == 0
  #@test winding_number(topo, P, SVector{3}(-1.0,-1.0,-1.0), cache3) |> round == 0
  @test winding_number(topo, P, SVector{3}(-1.0,-1.0,-1.0), cache.bvh) |> round == 0
  @test winding_number(topo, P, SVector{3}(-1.0,-1.0,-1.0)) |> round == 0
  @test winding_number(topo, P, SVector{3}(0.1,0.2,0.5), cache) |> round == 1 
  @test winding_number(topo, P, SVector{3}(0.1,0.2,0.5)) |> round == 1 

  @test HalfEdges.isorphan(floodfill(topo, P;verbose=true), 1) == false
  @test floodfill(topo, P) |> length == 1
  topo = Topology(vcat(polygons(topo), [[9, 10, 11]]))
  P = vcat(P, [SVector(0.5,3.0,-0.5), SVector(0.5,-1.0,3.5), SVector(0.5,-1.0,-1.0)])
  @test floodfill(topo, P) |> length == 3
end

@testset "load a mesh" begin
  topo, P = loadmesh("bunnylow.obj")
  @test typeof(topo) == Topology
  @test nfaces(topo) > 0
end

@testset ".Collision" begin
  Point = SVector{3,Float64}
  o,x,y,z = Point(0.,0,0), Point(1.,0,0), Point(0.,1,0), Point(0.,0,1)
  a,b = HalfEdges.Collision.segment_between(((o, x, y), 
                                       (0.2x+0.2y-0.5z, 0.2x+0.2y+0.5z, 0.2x - 0.5y+0.5z)),
                                       [(2, (1, 2)), (1, (1, 2))])
  @test a[3] == b[3] == 0.0
  @test HalfEdges.Collision.triangle_edges(o, x, y, Point(0.2, 0.2, 1.0),
                                              Point(0.2, 0.2, -1.0),
                                              Point(0.1, 0.1, -1.0)) |> first == true
  @test HalfEdges.Collision.triangle_edges(o, x, y, Point(0.2, 0.2, -0.1),
                                              Point(0.2, 0.2, -1.0),
                                              Point(0.1, 0.1, -1.0)) |> first == false
end

@testset "DifferentialGeometry" begin
  topo, P = cube()

  P = HalfEdges.poisson(HalfEdges.Δ(topo, P), zero.(P), P, [1,2,3,4])
  # minimal surface on x-z plane
  @test mapreduce(p->p[2], +, P) == 0.0
  @test mapreduce(p->p[1], *, P[5:end]) > 0.0
  @test mapreduce(p->p[1], *, P[5:end]) < 1.0

end

@testset "IslandFinder" begin
  @test length(find_islands([[1,2],[3,4],[5,6],[7,8]])) == 4
  @test length(find_islands([[1,2],[3,4],[5,6],[6,1]])) == 2
  @test length(find_islands([[1,2],[3,4],[5,6],[6,1],[4,5]])) == 1
end

include("BoundingVolumeTrees_runtests.jl")

