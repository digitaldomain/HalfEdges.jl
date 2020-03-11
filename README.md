[![Build Status](https://travis-ci.com/mewertd2/HalfEdges.jl.svg?branch=master)](https://travis-ci.com/mewertd2/HalfEdges.jl)
[![codecov.io](https://codecov.io/github/mewertd2/HalfEdges.jl/coverage.svg?branch=master)](https://codecov.io/github/mewertd2/HalfEdges.jl?branch=master)

# HalfEdges

The `HalfEdges` [Julia](http://julialang.org) package defines the `HalfEdges` type and 
methods/types to implement operations on the 
[halfedge data structure](https://en.wikipedia.org/wiki/Doubly_connected_edge_list).

At this time HalfEdges only supports immutable operations (read-only).  i.e. you can't currently modify a Topology once created.

## Creation

To represent the topology of a discrete polygon mesh, instance the `Topology` type.

    julia> topo = Topology([[1,2,3],[1,3,4],[1,4,5,2]],5)

or

    julia> topo = Topology([1,2,3,1,3,4,1,4,2])

a flat list of indices is assumed to be a manifold set of connected triangles.


## Handles

There are unique handle types for each simplex element.  These are essentially uniquely typed integer handles for the different simplex parts you might want to reference.
`EdgeHandle`, `VertexHandle`, and `FaceHandle`.  

`HalfEdgeHandle` is the most common return type for methods in the `HalfEdges` package.
`EdgeHandle` represents a unique edge, whereas there are two distinct `HalfEdgeHandle` per edge.

    julia> using HalfEdges: vertices, edges, face, halfedge, head, tail    

    julia> vertices(topo) |> first |> typeof
    HalfEdgeHandle

    julia> vertices(topo) |> first |> heᵢ->tail(topo, heᵢ) |> typeof
    VertexHandle

## HalfEdge conventions

A given `HalfEdge` points from `tail` to `head` `vertex` and is associated with the `face` on it's left.
A `Polygon` is considered to be "wound" in the counterclockwise direction.
`FaceHandle` and `VertexHandle` integer values will correspond to the polygon and vertex indices you created the `Topology` with, so long as you haven't called any methods that do any reindexing.  

    julia> using HalfEdges: halfedges, vertices, vertex, edges, edge, halfedge, face, head, tail, opposite, next

        2
       /|\
      / | \ 
     /  |  \
    3---1---4

    julia> topo = Topology([[1,2,3],[1,4,2]])

    # find the only halfedge pointing from vertex 1 to 2
    julia> h = Iterators.filter(h->head(topo, h)==VertexHandle(2), OneRing(topo, VertexHandle(1))) |> first
    3

    julia> tail(topo, h)
    1

    julia> head(topo, h)
    2

    julia> vertex(topo, h) == ans
    true

    julia> map( hᵢ->vertex(topo, hᵢ), Polygon(topo, face(topo, h)))
    3-element Array{VertexHandle,1}:
     3
     1
     2

    julia> h = opposite(topo, h)
    4

    julia> map( hᵢ->vertex(topo, hᵢ), Polygon(topo, face(topo, h)))
    3-element Array{VertexHandle,1}:
     1
     4
     2
 
## Iterating over components

A number of iterators are provided.  Unless stated, the iterates are `HalfEdgeHandles`

* `OneRing` iterates around the edge "spokes" radiating out from a vertex.  This iteration is done clockwise for performance reasons.
* `Rim` iterates around the "tire" edges of the polygons attached to a vertex.  This iteration is performned in a clockwise manner.
* `Polygon` iterates around a polgon in a counterclockwise manner.
* `BoundaryLoop` iterates around the boundary connected to the given halfedge.

Examples:

        2====6
       /+\  ||
      / + \ ||
     /  +  \||
    3+++1+++4
     \     /
      \   /
       \ /
        5 

    OneRing₁ => '+'
    Rim₁ => '-'

    julia> topo = Topology([[1,2,3],[1,4,2],[1,3,5,4],[2,4,6]])  # note we have a quad mixed in 

    julia> edges(topo, OneRing(topo, VertexHandle(1)))
    3-element Array{Tuple{VertexHandle,VertexHandle},1}:
     (1, 3)
     (1, 2)
     (1, 4)

    julia> edges(topo, Rim(topo, VertexHandle(1)))
    4-element Array{Tuple{VertexHandle,VertexHandle},1}:
     (2, 3)
     (4, 2)
     (3, 5)
     (5, 4)

    julia> BoundaryLoop(topo, first(OneRing(topo, VertexHandle(1))))

    julia> edges(topo, BoundaryLoop(topo, first(OneRing(topo, VertexHandle(2)))))
    5-element Array{Tuple{VertexHandle,VertexHandle},1}:
     (2, 6)
     (6, 4)
     (4, 5)
     (5, 3)
     (3, 2)

    julia> edges(topo, Polygon(topo, FaceHandle(1)))
    3-element Array{Tuple{VertexHandle,VertexHandle},1}:
     (3, 1)
     (1, 2)
     (2, 3)

## Navigating a mesh

Methods to move around by following `HalfEdgeHandle` are the most flexible ways to navigate.

Lets move from an edge on vertex 1 to vertex 6.

        2---6
       /^+  ^
      / + + +
     /  +  >+
    3++>1---4

We start on the halfedge 3->1, which is on the interior of the face (1,2,3).  
Use `next` to move ccw onto 1->2
`opposite` will jump across the edge to the halfedge 2->1 on face (1,4,2)
Now we pivot around vertex 2 over edge (2,4) using `next` and land on the halfedge 2->4 on face (4,6,2)
We finally move ccw  one more time with `next`
Finally check we are pointing at vertex 6

    # do it the long way
    julia> head(topo, next(topo, opposite(topo, next(topo, next(topo, opposite(topo, next(topo, HalfEdgeHandle(1))))))))
    6

    # use 1-arity methods to avoid repeating topo parameter excessively

    # make a function to pivot around a vertex over an edge
    julia> pivot = opposite(next(next))
    #62 (generic function with 1 method)

    # execute our plan to get from 3->1 to vertex 6
    julia> head(next)(topo, pivot(topo, opposite(next)(topo, HalfEdgeHandle(1))))
    6
    
    # use julia's function chaining macro to order operations in a more natural way.
    julia> (halfedge |> next |> opposite |> next |> next |> opposite |> next |> head)(topo, (3,1)) == ans
    true

## Collision Detection

There is support for accelerating collision queries using a Bounding Volume Hierarchy.

    # create a Collider first to cache the BVH
    julia> col = Collider(topo, P)

    # can query all the faces that overlap a bounding box
    julia> query_aabb(col, HalfEdges.BVH.AABB(SVector{3}(0.,0.,0.5), SVector{3}(1.0,1.0,1.0)))  

    # we need to provide our own method to actually collide the candidate faces
    julia> traingle_vs_triangle(a1,b1,c1, a2,b2,c2) = (true, "some collision data")  # provide your own

    # now we can do a self-intersection test on the entire mesh
    julia> hits = collide_self(col, triangle_vs_triangle)


## Project Information

### Contributing

Please read [CONTRIBUTING.md](./CONTRIBUTING.md) for details.

### Authors

* **Michael Alexander Ewert** - Developer - [Digital Domain](https://digitaldomain.com)

### License

This project is licensed under a modified Apache 2.0 license - see the [LICENSE](./LICENSE) file for details
