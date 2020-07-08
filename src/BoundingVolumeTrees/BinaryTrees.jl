#!me TODO:  
#===
#
#  Trait trick is probably needed to decomplect solve my type issues
f(a, storage::Contiguous) = ...
f(a, storage::Strided) = ...
f(a) = f(a, getstorage(a))

...

getstorage(::MyType) = Contiguous()
#
===#
#==
# iterators and path Reverse need to be sorted out
==#

export
EmptyNode,
Node,
Branch,
Leaves,
BinaryTree,
Traverse,
DepthFirst,
BreadthFirst,
BinarySearch,
Path,
ReversePath,
empty_node,
isempty_node,
left,
right,
data,
key,
value,
insert,
search,
tree,
leaves,
isleaf,
rotate_left,
rotate_right,
path,
reverse_path,
replace_node,
delete,
update


using Base.Iterators
using Base.Iterators:Reverse
using Base.Iterators:Drop

∂ = (f,x)-> (y...)->f(x,y...)

abstract type BinaryTree{T} end
const EmptyNode = Nothing
const empty_node = nothing

const Branch{T} = Union{BinaryTree{T},EmptyNode}

abstract type WrappedData{T} end

left(tree::EmptyNode) = nothing
right(tree::EmptyNode) = nothing
key(tree::EmptyNode) = nothing
value(tree::EmptyNode) = nothing
data(tree::EmptyNode) = nothing
keyval(tree::EmptyNode) = nothing

struct Node{T} <: BinaryTree{T}
  data::T
  left::Branch{T}
  right::Branch{T}
end

Node(data::T) where T = Node(data,empty_node,empty_node)

left(tree::Node{T}) where T = tree.left
right(tree::Node{T}) where T = tree.right
data(tree::Node{T}) where T = tree.data
key(tree::Node{T}) where T = key(tree.data)
value(tree::Node{T}) where T = value(tree.data)
keyval(tree::Node{T}) where T = tree.data
leaf(data::T) where T = Node(data)

key(d::Any) = d
value(d::Any) = d

internal_nodedata(::Type{T}) where T = one(T) 
internal_nodedata(::Type{T}) where {T<:AbstractArray} = T()
internal_nodedata(::Type{W}) where {T,W<:WrappedData{T}} = internal_nodedata(T) 

#  Traits  #

# a full tree always has either 0 or 2 children per branch
abstract type Fullness  end
struct FullTree <: Fullness end
struct NotFullTree <: Fullness end

fullness(::Type{Node{T}}) where {T} = NotFullTree
fullness(::Type{Node{W}}) where {T,W<:WrappedData{T}} = fullness(Node{T})

# AbstractTrees interface
import AbstractTrees.children
import AbstractTrees.printnode
import AbstractTrees.Leaves

children( n::Nothing ) = ()
children( n::BinaryTree{T} ) where T = (right(n),left(n))
printnode(io::IO, n::BinaryTree{T}) where T = show(data(n))
printnode(io::IO, n::Nothing ) = nothing


struct DepthFirst{T}
  tree::BinaryTree{T}
end

rest(a::Vector{T}) where T = Vector{T}(a[2:end])
rest(a::NTuple{T}) where T = a[2:end]

isempty_node(n::Nothing) = true
isempty_node(n) = false

isleaf(a::BinaryTree{T}) where T = (isempty_node(left(a)) && isempty_node(right(a)))

#==
abstract type Path end

struct ImmutablePath{Node{T}} <: Path
  v::Vector{PathElement{T}}
end

struct IndexPath{} end
==#

const PathElement{T} = Tuple{BinaryTree{T},Symbol}
const Path{T} = Vector{PathElement{T}}
#const ReversePath{T} = Reverse{Path{T}}
const ReversePath{T} = Union{Reverse{Path{T}},Drop{Reverse{Path{T}}}}  # awkward.

#==== General Traversal ===#
# Transducers should give better performance than Channel

function traverse(node::BinaryTree{T}, c::Channel, pre, post, inorder) where T
  if pre != nothing
    push!(c, pre(node))
  end

  traverse(left(node), c, pre, post, inorder)
  if inorder != nothing
    push!(c, inorder(node))
  end

  traverse(right(node), c, pre, post, inorder)
  if post != nothing
    push!(c, post(node))
  end
end

function traverse(node::Nothing, c::Channel, f...)
end

Traverse(tree::BinaryTree{T}, pre, post, inorder) where T = Channel() do c
  traverse(tree, c, pre, post, inorder)
end


#==== DepthFirst Traversal ===#
function Base.iterate(iter::DepthFirst{T}, children::Vector{BinaryTree{T}}) where T
  if length(children) == 0
    return nothing
  end
  child = first(children)
  (child, vcat(Vector{BinaryTree{T}}(filter(x->x!=nothing,[left(child), right(child)])), rest(children)))
end

function Base.iterate(iter::DepthFirst{T}) where T
  (iter.tree,Vector{BinaryTree{T}}(filter(x->x!=nothing,[left(iter.tree), right(iter.tree)])))
end

Base.length(iter::DepthFirst{T}) where T = reduce( (acc,_)->acc+1, iter; init=0)

Base.iterate(tree::BinaryTree{T}) where T = iterate(DepthFirst(tree))
Base.iterate(iter::BinaryTree{T}, children::Vector{BinaryTree{T}}) where T = iterate(DepthFirst(iter), children)
Base.length(iter::BinaryTree{T}) where T = length(DepthFirst{T}(iter))

#===== BreadthFirst Traversal ===#
struct BreadthFirst{T}
  tree::BinaryTree{T}
end

function Base.iterate(iter::BreadthFirst{T}, children::Vector{BinaryTree{T}}) where T
  if length(children) == 0
    return nothing
  end
  child = last(children)
  (child, vcat(Vector{BinaryTree{T}}(filter(x->x!=nothing,[right(child), left(child)])), children[1:end-1]))
end

function Base.iterate(iter::BreadthFirst{T}) where T
 (iter.tree,Vector{BinaryTree{T}}(filter(x->x!=nothing,[right(iter.tree), left(iter.tree)])))
end

Base.length(iter::BreadthFirst{T}) where T = length(DepthFirst(iter.tree))

#===== BinarySearch Traversal ===#
struct BinarySearch{T}
  tree::BinaryTree{T}
  comparison::Function
end

function pathfinder( k::T ) where T
  st->begin
    if key(k) == key(st)
      :found
    elseif key(k) < key(st)
      :left
    else
      :right
    end
  end
end

BinarySearch(tree::BinaryTree{T}, k::T) where {T} =
  BinarySearch{T}( tree, pathfinder(k) )

function Base.iterate(iter::BinarySearch{T}, subtree::Branch{T}) where T
  if subtree == empty_node
    return nothing
  end

  if iter.comparison(subtree) == :found
    ((subtree,:found),nothing)
  elseif iter.comparison(subtree) == :left
    ((subtree,:left), left(subtree))
  else
    ((subtree,:right), right(subtree))
  end
end

function Base.iterate(iter::BinarySearch{T}) where T
  iterate(iter, iter.tree)
end

Base.length(iter::BinarySearch{T}) where T = reduce( (acc,_)->acc+1, iter; init=0)

#height(t::BinaryTree{T}) where T = reduce(max,length.(∂(path,t).(key.(leaves(t)))))
height(::EmptyNode, h::Int) = h
height(t::Branch{T}, h::Int) where T = max(height(left(t),h+1), height(right(t), h+1))
height(t::BinaryTree{T}) where T = height(t, 0)

#===== Leaf Traversal ===#
Leaves( a::BinaryTree{T}, descend::Function = x->true ) where {T} = Channel(ctype=BinaryTree{T}) do c
  nodestack = Vector{BinaryTree{T}}()

  if descend(a)
    push!(nodestack,a)
  end

  while( !isempty(nodestack) )
    ancestor = pop!(nodestack)

    if( isleaf(ancestor) )
      push!(c, ancestor)
    else
      childL = left(ancestor)
      childR = right(ancestor)

      if childR != empty_node && descend(childR)
        push!(nodestack,childR)
      end
      if childL != empty_node && descend(childL)
        push!(nodestack,childL)
      end
    end
  end
end

leaves( a::BinaryTree{T}, descend::Function = x->true ) where T = Leaves( a, descend )

"""
  search nodes until ftest passes. returns the subtree rooted there
"""
search( ftest::Function, tree::BinaryTree{T} ) where T = Iterators.filter(ftest,tree) |> first

"""
  search nodes until we find the node with given value. returns the subtree rooted there
"""
search( k::T, tree::BinaryTree{T} ) where T = search(x->key(x)==k, tree)

"""
  search nodes until ftest passes, return the path to that node
"""
path( k::T, tree::BinaryTree{T} ) where T = path(BinarySearch(tree, k))
path( tree::BinaryTree{T}, k::T ) where T = path(k,tree)
path( bs::BinarySearch{T} ) where T =  (bs|>∂(collect,PathElement{T}))

"""
  reverse the path from top down to bottom up
"""
reverse( p::Path{T} ) where T = Reverse(p)
reverse( bs::BinarySearch{T} ) where T = Reverse(path(bs))

reverse_path = reverse∘path

"""
  returns a function to compare one node to another. for binary search. i.e. return :left if first < second
"""
branch_compare( d::S ) where {S} = st-> d < key(st) ? :left : :right

"""
  return a new tree with new node inserted on branch where data < parent.data decends to left
"""
insert( ::Type{NotFullTree}, t::N, d::T ) where {T, N<:BinaryTree{T}} =
  rbuild(reverse_path(BinarySearch(t, 
                                   branch_compare( key(d) ))),
         N(d,empty_node,empty_node))

function insert( ::Type{FullTree}, t::N, d::T ) where {T, N<:Node{T}} 
  rpath = reverse_path(BinarySearch(t, 
                                    branch_compare( key(d) )))

  # no entities are stored on branches, so create new Branch with existing Leaf right and new Leaf left
  cleaf = (take(rpath,1) |> first)[1]
  rpath = drop(rpath,1)

  twig = rebuild_leaf( FullTree, :both, cleaf, d )
  rbuild( rpath, twig )
end

insert( t::N, d::T ) where {T, N<:BinaryTree{T}} = insert( fullness(N), t, d )
#insert( t::N, d::T ) where {T, N<:Node{T}} = insert( fullness(N), t, d )
#insert( t::N, d::T ) where {T, N<:BinaryTree{T}} = insert( fullness(N), t, d )
"""
  rebuild a tree given a path.  path will be created, all other branches preserved
"""
build( path::Path{T} ) where T =  rbuild( Reverse(path) )

rebuild_leaf( ::Type{NotFullTree}, direction::Symbol, parent::N, d::T ) where {N,T} =
  rebuild_node( Val(direction), empty_node, empty_node, parent, d )

function rebuild_leaf( ::Type{FullTree}, direction::Symbol, fatherbrother::N, d::T ) where {N,T}
  newleaf = rebuild_node( Val(direction), empty_node, empty_node, fatherbrother, d )
  twigd = rebuild_parent_data( Val(:both), keyval(newleaf), keyval(fatherbrother), nothing ) 
  # well now, this is incestuous. the last fatherbrother is ignored usually.
  #!me maybe that last fatherbrother should be nothing?  then change rebuild_node to accept Branch there
  #!me might result in ambiguous dispatch for wrapped type?
  twig = rebuild_node(Val(:both), newleaf, fatherbrother, fatherbrother, twigd)
end

#WIP
rebuild_parent_data(dir,#::Union{Val{:both}, Val{:left}, Val{:right}, Val{:either}}, 
                    lcd::Union{T,Nothing}, rcd::Union{T,Nothing}, 
                    parent_d::T ) where T = parent_d

rebuild_parent_data(dir, #::Union{Val{:both}, Val{:left}, Val{:right}, Val{:either}}, 
                    lcd::Union{T,Nothing}, rcd::Union{T,Nothing}, 
                    _::Nothing ) where {T<:Union{AbstractString,Number,AbstractArray}} = internal_nodedata(T)

rebuild_node(::Union{Val{:both}, Val{:left}, Val{:right}}, 
             lc::Branch{T}, rc::Branch{T}, 
             n::Node{T}, d::T ) where T = 
  Node(d, lc, rc) 

function rebuild_node( direction::Symbol, lc::Branch{T}, rc::Branch{T}, n::Node{T} ) where T
  d = rebuild_parent_data( Val(direction), keyval(lc), keyval(rc), keyval(n) ) 
  rebuild_node( Val(direction), lc, rc, n, d )
end

function rebuild_node(direction::Symbol, child::Branch{T}, n::BinaryTree{T} ) where T
  if direction == :left
    rebuild_node( direction, child, right(n), n )
  else
    rebuild_node( direction, left(n), child, n )
  end
end

update(n::Branch{T}) where T = n
#WIP

"""
  rebuild a tree given a reversed (bottom up) path.  path will be created, all other branches preserved
  to implement a new type of BinaryTree, add a new method with this signature for rebuilding
"""
function rbuild( rpath::ReversePath{T}, leaf::Branch{T} ) where T
  reduce((st,(n, direction)) -> rebuild_node( direction, st, n ),
         rpath ;
         init=update(leaf)) 
end

rbuild( rpath::ReversePath{T} ) where T = rbuild(drop(rpath,1), (first∘first)(rpath))

"""
  replace the deepest node in dest path with source subtree
"""
function replace_node( dest::ReversePath{T}, source::Branch{T} ) where T
  rpath = drop(dest, 1)
  rbuild(rpath,source)
end

function delete( rpath::ReversePath{T} ) where T
  node,dir = first(rpath)
  dir == :found || return first(last(rpath)) 
  lc = left(node)
  rc = right(node)

  if lc == empty_node && rc == empty_node
    replace_node( rpath, empty_node )
  elseif lc == empty_node
    replace_node(rpath, rc)
  elseif rc == empty_node
    replace_node(rpath, lc)
  else  # both children present
    # choose left or right?  right for now ( right then take leftmost )
    # the new root of this subtree will be leftmost leaf of right branch
    goleft = x->:left
    successor_path = reverse_path(BinarySearch(rc,goleft))
    successor = (first∘first)(successor_path)
    #rtree = replace_node(successor_path, empty_node)
    rtree = replace_node(successor_path, right(successor))
    #!me NG
    #replace_node(rpath, Node(data(predecessor),lc,rtree))
    #newsource = rebuild_node(:left, lc,  
    newelder = rebuild_node(:both, lc, rtree, successor)
    replace_node(rpath, newelder)
  end

end


delete( t::BinaryTree{T}, key::T ) where T = delete(reverse_path(key,t))
delete( key::T, t::BinaryTree{T} ) where T = delete(t, key)
delete( t::BinaryTree{W}, d::T ) where {T,W<:WrappedData{T}} = delete(t, W(d))

tree( d::T ) where T = Node(d,empty_node,empty_node)

tree( dc::Vector{T} ) where T = reduce(insert,rest(dc); init = tree(first(dc)))

function rotate_right( parent::Node{T} ) where T
  child = left(parent)
  Node(data(child), left(child), Node(data(parent), right(child), right(parent)))
end

function rotate_left( parent::Node{T} ) where T
  child = right(parent)
  Node(data(child), Node(data(parent), left(parent), left(child)), right(child))
end

update( t::BinaryTree{T}, data::T, newdata::T ) where T = insert(delete(t, key(data)), newdata)

include("BalancedTrees.jl")

