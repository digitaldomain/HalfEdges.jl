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
end

@testset "3x3grid" begin
  cell(i,j) = map(v->v+(j-1)*4, [i,i+4,i+1,i+1,i+4,i+5])
  topo = Topology(Iterators.flatten([cell(1,1), cell(2,1), cell(3,1), 
                                     cell(1,2), cell(2,2), cell(3,2), 
                                     cell(1,3), cell(2,3), cell(3,3)])|>collect)
  @test sort(he_.boundary_interior(topo)) == [6,7,10,11]

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


