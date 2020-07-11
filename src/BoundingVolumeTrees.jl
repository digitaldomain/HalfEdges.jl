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

#!me this one is optimized, but lacking features like balancing
#!me also the design is a bit grim
# however, it works, so it's the one we use in HalfEdges for now
include("BoundingVolumeTrees/DynamicAABBTrees.jl")

end
