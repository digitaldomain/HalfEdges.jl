module IndexedBinaryTree

struct IndexedNode{T}
  data::T
  left::Int
  right::Int
end

const EmptyIndexedNode = -1

struct NodeVector{T}
  nodes::Vector{IndexedNode{T}} 
  root::Int
end

IndexedNode(a::T) where T = IndexedNode{T}(a, EmptyIndexedNode, EmptyIndexedNode)

struct IndexedBinaryTree{T} <: BinaryTree{T}
  nodes::NodeVector{T}
  node::Int
end

Base.size(A::NodeVector{T}) where T = size(A.nodes)
Base.getindex(A::NodeVector{T}, i::Number) where {T} = A.nodes[i]
Base.firstindex(A::NodeVector{T}) where T = A.root

itree(ibt::IndexedBinaryTree{T}, node::Int) where T = IndexedBinaryTree{T}(ibt.nodes, node)
itree(a::T) where T = IndexedBinaryTree{T}(NodeVector{T}([IndexedNode(a)], 1), 1)

left(tree::IndexedBinaryTree{T}) where T = isempty_node(tree.nodes[tree.node].left) ? nothing : itree(tree, tree.nodes[tree.node].left)
right(tree::IndexedBinaryTree{T}) where T = isempty_node(tree.nodes[tree.node].right) ? nothing : itree(tree, tree.nodes[tree.node].right)


isempty_node(tree::IndexedBinaryTree{T}) where T = tree.node == EmptyIndexedNode
isempty_node(i::Int) = i == EmptyIndexedNode

data(tree::IndexedBinaryTree{T}) where T = tree.nodes[tree.node].data

end
