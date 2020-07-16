module BoundingVolumeTrees

using StaticArrays
import StaticArrays.insert
const Point{T} = SVector{3,T}
const Vector3 = Point

# some WIP tree structures intended to replace the older DynamicAABBTree
include("BoundingVolumeTrees/BinaryTrees.jl")
include("BoundingVolumeTrees/IndexedBinaryTrees.jl") 
include("BoundingVolumeTrees/BalancedTrees.jl")
include("BoundingVolumeTrees/AABBs.jl")
include("BoundingVolumeTrees/AABBTrees.jl")


#!me ideally this isn't needed.  attempt to have a pool of leaves at indices 1:nleaves reserved
function index_entities(aabbtree::B) where {T,S,B<:Union{IndexedBinaryTree{AVLData{AABBNodeData{T,S}}},
                                                         IndexedBinaryTree{AABBNodeData{T,S}}}}
  [x[2] for x in sort(map(n->(entity(value(data(n))),index(n)),Leaves(bt)), by=first)]
end



#!me this one is optimized, but lacking features like balancing
#!me also the design is a bit grim
# however, it works, so it's the one we use in HalfEdges for now
include("BoundingVolumeTrees/DynamicAABBTrees.jl")

end
