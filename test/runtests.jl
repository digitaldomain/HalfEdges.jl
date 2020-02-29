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

@testset "one tri" begin
  topo = HalfEdges.Topology([[1,2,3]], 3)
  @test 1 == 1  # just make sure it didn't fail out
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
end


@testset "unbroken ring" begin
  topo = Topology([[1,2,3],[1,3,4],[1,4,2]],4)

  @test length(collect(HalfEdges.OneRing(topo,VertexHandle(1)))) == 3
  @test length(collect(HalfEdges.OneRing(topo,VertexHandle(4)))) == 3
  @test reduce( (acc,iv)->acc && iv == VertexHandle(1), [HalfEdges.tail(topo,ih) for ih in HalfEdges.OneRing(topo,VertexHandle(1))][1:end], init=true)
  @test nhalfedges(topo) == length(edges(topo))*2
end

@testset "box48" begin
  he_ = HalfEdges
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
