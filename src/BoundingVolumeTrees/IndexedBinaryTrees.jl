export 
IndexedBinaryTree,
itree,
index,
ileaves,
itraversewhen

const EmptyIndexedNode = -1

struct IndexedNode{T}
  left::Int
  right::Int
  next::Int # next can mean parent or next free node
  data::T
  IndexedNode{T}() where T = new{T}(EmptyIndexedNode, EmptyIndexedNode, EmptyIndexedNode)
  IndexedNode{T}(l,r) where T = new{T}(l, r, EmptyIndexedNode)
  IndexedNode{T}(l,r,n) where T = new{T}(l, r, n)
  IndexedNode{T}(l,r,n,d::T) where T = new{T}(l, r, n, d)
end

IndexedNode(a::T) where T = IndexedNode{T}(EmptyIndexedNode, EmptyIndexedNode, EmptyIndexedNode, a)
IndexedNode(l,r,n,d::T) where T = IndexedNode{T}(l, r, n, d)

mutable struct NodeVector{T}
  nodes::Vector{IndexedNode{T}} 
  root::Int
  nextfree::Int
  last::Int
  #nextfreeleaf::Int  # can optionally reserve first block of indexes for leaves.  Enables direct indexing
end

isleaf(n::IndexedNode{T}) where T = EmptyIndexedNode == n.left && 
                                    EmptyIndexedNode == n.right

parent(n::IndexedNode{T}) where T = n.next
setparent(n::IndexedNode{T}, iparent) where T = IndexedNode{T}(n.left, n.right, iparent, n.data)
setnext(n::IndexedNode{T}, inext) where T = setparent(n, inext)

Base.length(nv::NodeVector{T}) where T = length(nv.nodes)

struct IndexedBinaryTree{T} <: BinaryTree{T}
  nodes::NodeVector{T}
  node::Int
end

IndexedBinaryTree{T}(d::T, ::Nothing, ::Nothing) where T = itree(d)
IndexedBinaryTree(d::T) where T = itree(d)
IndexedBinaryTree{T}(d::T) where T = itree(d)

isroot(t::IndexedBinaryTree{T}) where T = t.nodes.root == t.node

function NodeVector{T}(sz::Int) where T
  nodes = [IndexedNode{T}(EmptyIndexedNode, EmptyIndexedNode, i) for i in vcat(2:sz, EmptyIndexedNode)]
  NodeVector{T}(nodes, EmptyIndexedNode, 1, sz)
end

index(tree::IndexedBinaryTree{T}) where T = tree.node
node(tree::IndexedBinaryTree{T}) where T = tree.nodes[tree.node]

root(tree::IndexedBinaryTree{T}) where T = IndexedBinaryTree{T}(tree.nodes, tree.nodes.root)

Base.getindex(A::IndexedBinaryTree{T}, i::Number) where {T} = IndexedBinaryTree{T}(A.nodes, i)

function ileaves(t::IndexedBinaryTree{T}) where T
  L = Vector{IndexedNode{T}}(undef, (length(t.nodes)>>1)+1)
  resize!(L, 0);
  ileaves(t.nodes, t.node, L);
  return L
end

function ileaves( n::NodeVector{T}, i::Int, L::Vector{IndexedNode{T}} ) where T
  nc = 0
  if n.nodes[i].left != EmptyIndexedNode
    nc = 1
    ileaves(n, n.nodes[i].left, L)
  end
  if n.nodes[i].right != EmptyIndexedNode
    nc += 1
    ileaves(n, n.nodes[i].right, L)
  end

  if nc == 0
    push!(L, n.nodes[i])
  end
end


function itraversewhen(when::Function, node::IndexedBinaryTree{T}, L::Vector{R}, 
                       pre::Union{Function, Nothing}, 
                       post::Union{Function, Nothing}, 
                       inorder::Union{Function, Nothing}) where {T, R}
  if !when(node)
    return
  end

  if pre != nothing
    push!(L, pre(node))
  end

  node.nodes[node.node].left != EmptyIndexedNode && itraversewhen(when, left(node), L, pre, post, inorder)

  if inorder != nothing
    push!(L, inorder(node))
  end

  node.nodes[node.node].right != EmptyIndexedNode && itraversewhen(when, right(node), L, pre, post, inorder)
  if post != nothing
    push!(L, post(node))
  end
end

function itraversewhen(when::Function, t::IndexedBinaryTree{T}, returnType::Type, 
                       pre, post, inorder = nothing) where T
  L = Vector{returnType}(undef, length(t.nodes))
  resize!(L, 0);
  itraversewhen(when, t, L, pre, post, inorder);
  return L
end

function itraversewhen(when::Function, node::IndexedBinaryTree{T}, L::Vector{R}, 
                       pre::Function) where {T, R}
  if !when(node)
    return
  end

  push!(L, pre(node))

  node.nodes[node.node].left != EmptyIndexedNode && itraversewhen(when, left(node), L, pre)
  node.nodes[node.node].right != EmptyIndexedNode && itraversewhen(when, right(node), L, pre)
end

function itraversewhen(when::Function, t::IndexedBinaryTree{T}, returnType::Type, 
                       pre::Function ) where T
  L = Vector{returnType}(undef, length(t.nodes))
  resize!(L, 0);
  itraversewhen(when, t, L, pre);
  return L
end

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

import Base.resize!
function resize!( nv::NodeVector{T}, sz::Int ) where T
  oldsz = length(nv.nodes)
  if sz < oldsz
    return nv
  end

  if oldsz*2 > sz
    sz = oldsz*2
  end
  resize!(nv.nodes, sz)
  newnodes = [IndexedNode{T}(EmptyIndexedNode, EmptyIndexedNode, i) 
              for i in vcat((oldsz+2):sz, EmptyIndexedNode)]
  nv.nodes[oldsz+1:sz] = newnodes
  lastn = nv.nodes[nv.last]
  nv.nodes[nv.last] = IndexedNode{T}(lastn.left, lastn.right, oldsz+1, lastn.data)

  if nv.nextfree == EmptyIndexedNode
    nv.nextfree = oldsz+1
  end

  nv.last = sz
  return nv
end

"""
    reserve_leaves!(nv, n) 

For internal use.
reserve the next n nodes for allocating to leave in order of insertion

this allows direct indexing of leaves with implicit in-order indices
uses any exising free nodes, so if those aren't ordered then this will not be true
"""
#==
function reserve_leaves!( nv::NodeVector{T}, n:Int ) where T
  nv.nextfreeleaf = nv.last
  oldsz = length(nv.nodes)
  resize!(nv, n+oldsz+1)
  nv.nextfree = n
  nv.nodes[nv.nextfreeleaf+n-1].next = EmptyIndexedNode
   
end
==#

Base.size(A::NodeVector{T}) where T = size(A.nodes)
Base.getindex(A::NodeVector{T}, i::Number) where {T} = A.nodes[i]
Base.setindex!(A::NodeVector{T}, n::IndexedNode{T}, i::Number) where {T} = A.nodes[i] = n
Base.firstindex(A::NodeVector{T}) where T = A.root

itree(ibt::IndexedBinaryTree{T}, node::Int) where T = IndexedBinaryTree{T}(ibt.nodes, node)
itree(a::T) where T = IndexedBinaryTree{T}(NodeVector{T}([IndexedNode(a)], 1, EmptyIndexedNode, 1), 1)
itree(dc::V) where {T, V<:AbstractVector{T}} = reduce(insert, Iterators.drop(dc,1); init = itree(first(dc)))

key(tree::IndexedBinaryTree{T}) where T = key(tree.nodes[tree.node].data)

left(tree::IndexedBinaryTree{T}) where T = isempty_node(tree.nodes[tree.node].left) ? nothing : itree(tree, tree.nodes[tree.node].left)

right(tree::IndexedBinaryTree{T}) where T = isempty_node(tree.nodes[tree.node].right) ? nothing : itree(tree, tree.nodes[tree.node].right)

parent(tree::IndexedBinaryTree{T}) where T = isempty_node(tree.nodes[tree.node].parent) ? nothing : itree(tree, parent(tree.nodes[tree.node]))

isempty_node(tree::IndexedBinaryTree{T}) where T = tree.node == EmptyIndexedNode
isempty_node(i::Int) = i == EmptyIndexedNode

data(tree::IndexedBinaryTree{T}) where T = tree.nodes[tree.node].data
keyval(tree::IndexedBinaryTree{T}) where T = data(tree)
data(n::IndexedNode{T}) where T = n.data
keyval(n::IndexedNode{T}) where T = n.data

function search(stree::IndexedBinaryTree{T}, k::T) where T
  hit = nothing
  while true
    if k == key(stree)
      hit = (stree, :found)
      break
    elseif k < key(stree)
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
  return hit
end

#==
function insert(::Type{NotFullTree}, t::IndexedBinaryTree{T}, d::T) where T
  nv = t.nodes
  if t.nodes.nextfree == EmptyIndexedNode
    nv = resize(nv, length(nv)+1)
  end
  inode = nv.nextfree
  nv_nextfree = nv[inode].next

  hit = search(t, d)
  hit = isnothing(hit) ? (t, :found) : hit

  iparent = first(hit).node
  p = nv[iparent]
  if hit[2] == :left
    p = IndexedNode{T}(inode, p.right, parent(p), p.data)
  else
    p = IndexedNode{T}(p.left, inode, parent(p), p.data)
  end
  
  cpnodes = copy(nv.nodes)
  cpnodes[inode] = IndexedNode{T}(EmptyIndexedNode, EmptyIndexedNode, iparent, d)
  cpnodes[iparent] =  p
  IndexedBinaryTree{T}(NodeVector{T}(cpnodes, nv.root, nv_nextfree, nv.last), nv.root)
end

import Base.insert!

insert!( t::N, d::T ) where {T, N<:BinaryTree{T}} = insert!( fullness(N), t, d )

function insert!(::Type{NotFullTree}, tree::IndexedBinaryTree{T}, d::T) where T
  nv = tree.nodes
  if tree.nodes.nextfree == EmptyIndexedNode
    nv = resize!(nv, length(nv)+1)
  end
  inode = nv.nextfree
  nv_nextfree = nv[inode].next

  hit = search(tree, d);
  hit = isnothing(hit) ? (tree, :found) : hit

  iparent = first(hit).node
  p = nv[iparent]
  if hit[2] == :left
    p = IndexedNode{T}(inode, p.right, parent(p), p.data)
  else
    p = IndexedNode{T}(p.left, inode, parent(p), p.data)
  end
  
  cpnodes = nv.nodes
  cpnodes[inode] = IndexedNode{T}(EmptyIndexedNode, EmptyIndexedNode, iparent, d)
  cpnodes[iparent] =  p
  IndexedBinaryTree{T}(NodeVector{T}(cpnodes, nv.root, nv_nextfree, nv.last), nv.root)
end
==# 

function rebuild_node!(dir::Union{Val{:both}, Val{:left}, Val{:right}}, 
             lc::Branch{T}, rc::Branch{T}, 
             n::IndexedBinaryTree{T}, d::T ) where T

  nv = n.nodes
  inode = n.node
  if isnothing(lc)
    ilc = EmptyIndexedNode
  elseif dir != Val{:right}() && isleaf(lc)
    # copy leaf into parents vector
    if nv.nextfree == EmptyIndexedNode
      nv = resize!(nv, length(nv)+2)  # adding two for case of 2 children
    end

    #==
    #!me possible orphaning? then we need something like this, but then that should looks like
    #!me else case after this
    #!me can happen for full trees (:both) or rotations we move a node from this tree to this tree?
    if !isroot(lc)
      nv[index(lc)] = setnext(node(lc), nv.nextfree)
      nv.nextfree = index(lc)
    end
    ==#

    ilc = nv.nextfree   # if not new index it will orphan current lc. else nv[index(lc)] = setparent(node(lc), nv.nextfree);  nv.nextfree = index(lc);  # maybe put in free() method
    nv.nextfree = nv[ilc].next
    nv[ilc] = setparent(node(lc), inode)
  else
    # child and parents vector should be the same
    ilc = index(lc)
    nv[ilc] = setparent(node(lc), inode)
  end

  if isnothing(rc)
    irc = EmptyIndexedNode
  elseif dir != Val{:left}() && isleaf(rc)
    # copy leaf into parents vector
    if nv.nextfree == EmptyIndexedNode
      nv = resize!(nv, length(nv)+1)
    end
    irc = nv.nextfree    
    nv.nextfree = nv[irc].next
    nv[irc] = setparent(node(rc), inode)
  else
    # child and parents vector should be the same
    irc = index(rc)
    nv[irc] = setparent(node(rc), inode)
  end

  imabastard = EmptyIndexedNode
  tn = IndexedNode{T}(ilc, irc, imabastard, d)
  nv[inode] = tn
  IndexedBinaryTree{T}(nv, inode)
end

rebuild_node(dir::Union{Val{:both}, Val{:left}, Val{:right}}, 
             lc::Branch{T}, rc::Branch{T}, 
             n::IndexedBinaryTree{T}, d::T ) where T = 
  rebuild_node!(dir, lc, rc, n, d)

