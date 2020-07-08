const EmptyIndexedNode = -1

struct IndexedNode{T}
  left::Int
  right::Int
  next::Int
  data::T
  IndexedNode{T}() where T = new{T}(EmptyIndexedNode, EmptyIndexedNode, EmptyIndexedNode)
  IndexedNode{T}(l,r) where T = new{T}(l, r, EmptyIndexedNode)
  IndexedNode{T}(l,r,n) where T = new{T}(l, r, n)
  IndexedNode{T}(l,r,n,d::T) where T = new{T}(l, r, n, d)
end

IndexedNode(a::T) where T = IndexedNode{T}(EmptyIndexedNode, EmptyIndexedNode, EmptyIndexedNode, a)
IndexedNode(l,r,n,d::T) where T = IndexedNode{T}(l, r, n, d)

#mutable struct NodeVector{T}
struct NodeVector{T}
  nodes::Vector{IndexedNode{T}} 
  root::Int
  nextfree::Int
  last::Int
end

Base.length(nv::NodeVector{T}) where T = length(nv.nodes)

struct IndexedBinaryTree{T} <: BinaryTree{T}
  nodes::NodeVector{T}
  node::Int
end

IndexedBinaryTree{T}(d::T, ::Nothing, ::Nothing)  where T = itree(d)

function NodeVector{T}(sz::Int) where T
  nodes = [IndexedNode{T}(EmptyIndexedNode, EmptyIndexedNode, i) for i in vcat(2:sz, EmptyIndexedNode)]
  NodeVector{T}(nodes, EmptyIndexedNode, 1, sz)
end

root(tree::IndexedBinaryTree{T}) where T = IndexedBinaryTree{T}(tree.nodes, tree.nodes.root)

function resize( nv::NodeVector{T}, sz::Int ) where T
  oldsz = length(nv.nodes)
  if sz < oldsz
    return nv
  end

  if oldsz*2 > sz
    sz = oldsz*2
  end
  newnodes = [IndexedNode{T}(EmptyIndexedNode, EmptyIndexedNode, i) 
              for i in vcat((oldsz+2):sz, EmptyIndexedNode)]
  nodes = vcat(nv.nodes, newnodes)
  lastn = nodes[nv.last]
  nodes[nv.last] = IndexedNode{T}(lastn.left, lastn.right, oldsz+1, lastn.data)
  NodeVector{T}(nodes, nv.root, oldsz+1, sz) 
end


Base.size(A::NodeVector{T}) where T = size(A.nodes)
Base.getindex(A::NodeVector{T}, i::Number) where {T} = A.nodes[i]
Base.setindex!(A::NodeVector{T}, n::IndexedNode{T}, i::Number) where {T} = A.nodes[i] = n
Base.firstindex(A::NodeVector{T}) where T = A.root

itree(ibt::IndexedBinaryTree{T}, node::Int) where T = IndexedBinaryTree{T}(ibt.nodes, node)
itree(a::T) where T = IndexedBinaryTree{T}(NodeVector{T}([IndexedNode(a)], 1, EmptyIndexedNode, 1), 1)

fullness(::Type{IndexedBinaryTree{T}}) where {T} = NotFullTree
fullness(::Type{IndexedNode{T}}) where {T} = NotFullTree
fullness(::Type{IndexedBinaryTree{W}}) where {T,W<:WrappedData{T}} = fullness(IndexedNode{T})

key(tree::IndexedBinaryTree{T}) where T = key(tree.nodes[tree.node].data)
left(tree::IndexedBinaryTree{T}) where T = isempty_node(tree.nodes[tree.node].left) ? nothing : itree(tree, tree.nodes[tree.node].left)
right(tree::IndexedBinaryTree{T}) where T = isempty_node(tree.nodes[tree.node].right) ? nothing : itree(tree, tree.nodes[tree.node].right)


isempty_node(tree::IndexedBinaryTree{T}) where T = tree.node == EmptyIndexedNode
isempty_node(i::Int) = i == EmptyIndexedNode

data(tree::IndexedBinaryTree{T}) where T = tree.nodes[tree.node].data

function insert(tree::IndexedBinaryTree{T}, d::T) where T
  nv = tree.nodes
  if tree.nodes.nextfree == EmptyIndexedNode
    nv = resize(nv, length(nv)+1)
  end
  inode = nv.nextfree
  nv_nextfree = nv[inode].next

  hit = BinarySearch(tree, d) |> collect |> last
  #==
  stree = tree
  hit = (stree, :found)
  while true
    if d == key(stree)
      hit = (stree, :found)
      break
    elseif d < key(stree)
      if isnothing(left(stree))
        hit = (stree, :left)
        break
      else
        stree = left(stree)
      end
    else 
      if isnothing(right(stree))
        hit = (stree, :right)
        break
      else
        stree = right(stree)
      end
    end
  end
  ==#

  iparent = first(hit).node
  parent = nv[iparent]
  if hit[2] == :left
    parent = IndexedNode{T}(inode, parent.right, parent.next, parent.data)
  else
    parent = IndexedNode{T}(parent.left, inode, parent.next, parent.data)
  end
  
  cpnodes = copy(nv.nodes)
  cpnodes[inode] = IndexedNode{T}(EmptyIndexedNode, EmptyIndexedNode, nv[inode].next, d)
  cpnodes[iparent] =  parent
  IndexedBinaryTree{T}(NodeVector{T}(cpnodes, nv.root, nv_nextfree, nv.last), nv.root)
end

import Base.import!
function insert!(tree::IndexedBinaryTree{T}, d::T) where T
  nv = tree.nodes
  if tree.nodes.nextfree == EmptyIndexedNode
    nv = resize(nv, length(nv)+1)
  end
  inode = nv.nextfree
  nv_nextfree = nv[inode].next

  #hit = BinarySearch(tree, d) |> collect |> last
  stree = tree
  hit = (stree, :found)
  while true
    if d == key(stree)
      hit = (stree, :found)
      break
    elseif d < key(stree)
      if isnothing(left(stree))
        hit = (stree, :left)
        break
      else
        stree = left(stree)
      end
    else 
      if isnothing(right(stree))
        hit = (stree, :right)
        break
      else
        stree = right(stree)
      end
    end
  end

  iparent = first(hit).node
  parent = nv[iparent]
  if hit[2] == :left
    parent = IndexedNode{T}(inode, parent.right, parent.next, parent.data)
  else
    parent = IndexedNode{T}(parent.left, inode, parent.next, parent.data)
  end
  
  cpnodes = nv.nodes
  cpnodes[inode] = IndexedNode{T}(EmptyIndexedNode, EmptyIndexedNode, nv[inode].next, d)
  cpnodes[iparent] =  parent
  IndexedBinaryTree{T}(NodeVector{T}(cpnodes, nv.root, nv_nextfree, nv.last), nv.root)
end


#==
rebuild_node(::Union{Val{:both}, Val{:left}, Val{:right}}, 
             lc::Nothing, rc::Nothing, 
             n::IndexedBinaryTree{T}, d::T ) where T =
  n.nodes

rebuild_node(::Union{Val{:both}, Val{:left}, Val{:right}}, 
             lc::IndexedBinaryTree{T}, rc::Nothing, 
             n::IndexedBinaryTree{T}, d::T ) where T =
  resize( 

rebuild_node(::Union{Val{:both}, Val{:left}, Val{:right}}, 
             lc::Nothing, rc::Nothing, 
             n::IndexedBinaryTree{T}, d::T ) where T =
  itree(d) 

rebuild_node(::Union{Val{:both}, Val{:left}, Val{:right}}, 
             lc::Nothing, rc::Nothing, 
             n::IndexedBinaryTree{T}, d::T ) where T =
  itree(d) 
  ==#
