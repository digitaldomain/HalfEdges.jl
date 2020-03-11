export
AVLData,
BalancedTree,
BalanceData,
AVLTree,
avltree,
avldata,
weigh,
height,
build_wrapped_node

abstract type BalanceTrait end
struct Unbalanced <: BalanceTrait end
struct AVLBalanced <: BalanceTrait end
struct WAVLBalanced <: BalanceTrait end
struct RAVLBalanced <: BalanceTrait end

abstract type BalanceData{T} <: WrappedData{T} end
struct AVLData{T} <: BalanceData{T}
  kv::T
  rank::Int32
end

wrapped_data(n::AVLData{T}) where T = n.kv

build_wrapped_node(::Type{W},d::T) where {W<:WrappedData,T} = Node(W(d), empty_node,empty_node)
build_wrapped_node(::Type{T},d::T) where {T} = Node(d, empty_node,empty_node)
build_wrapped_node(::Type{W},d::T, l, r) where {W<:WrappedData,T} = Node(W(d), l, r)
build_wrapped_node(::Type{T},d::T, l, r) where {T} = Node(d, l, r)

AVLData(kv::T, h::S) where {T,S<:Integer} = AVLData(kv,Int32(h))
AVLData(kv::T) where T = AVLData(kv,Int32(0))
AVLData{T}(kv::T) where T = AVLData(kv,Int32(0))

const AVLNode{T} = Node{AVLData{T}}

balance_trait(::Type) = Unbalanced()
balance_trait(::Type{EmptyNode}) = Unbalanced()
balance_trait(::Type{<:AVLNode{T}}) where T = AVLBalanced()

avltree( d::T ) where T = tree( AVLData(d) )
avltree( dc::Vector{T} ) where T = tree( AVLData.(dc) )
avldata( kv::T ) where T = AVLData(kv)

function rotate_right( c::AVLNode{T} ) where T
  b = left(c)
  indata = rebuild_parent_data( Val(:either), keyval(right(b)), keyval(right(c)), keyval(c) ) 
  cdata = AVLData(indata, height(right(b),right(c)))
  c = Node(cdata,right(b),right(c))

  indata = rebuild_parent_data( Val(:either), keyval(left(b)), keyval(c), keyval(b) ) 
  Node(AVLData(indata, height(left(b),c)), left(b), c)
end

function rotate_left( c::AVLNode{T} ) where T
  b = right(c)
  indata = rebuild_parent_data( Val(:either), keyval(left(c)), keyval(left(b)), keyval(c) ) 
  c = Node(AVLData(indata, height(left(c),left(b))),left(c),left(b))

  indata = rebuild_parent_data( Val(:either), keyval(c), keyval(right(b)), keyval(b) ) 
  Node(AVLData(indata, height(c,right(b))), c, right(b))
end

balance( n::T ) where T = balance( n, balance_trait(T) ) 

balance( n::T, ::Unbalanced ) where T = n

function balance( st::AVLNode{T}, trait::AVLBalanced ) where T
  scale = weigh(st)
  if scale > 1
    # left heavy
    subscale = weigh(left(st))
    if subscale > 0
      # left child heavy
      rotate_right(st)
    else
      # right child heavy
      rotate_right(Node(data(st), rotate_left(left(st)), right(st)))
    end
  elseif scale < -1
    # right heavy
    subscale = weigh(left(st))
    if subscale > 0
      # left child heavy
      rotate_left(Node(data(st), left(st), rotate_right(right(st))))
    else
      # right child heavy
      rotate_left(st)
    end
  else
    st
  end
end

function rebuild_node(::Val{:left}, 
                      child::Branch{AVLData{T}}, 
                      _::Branch{AVLData{T}}, 
                      n::Node{AVLData{T}}, d::T ) where T
  avldata = AVLData{T}(d, max(height(right(n)),height(child))+1)
  balance(Node(avldata, child, right(n)))
end

function rebuild_node(::Val{:right}, 
                      _::Branch{AVLData{T}}, 
                      child::Branch{AVLData{T}}, 
                      n::Node{AVLData{T}}, d::T ) where T
  avldata = AVLData{T}(d, max(height(left(n)),height(child))+1)
  balance(Node(avldata,left(n),child))
end

function rebuild_node(::Val{:both}, 
                      lc::Branch{AVLData{T}}, 
                      rc::Branch{AVLData{T}}, 
                      n::Node{AVLData{T}}, d::T ) where T
  avldata = AVLData{T}(d, max(height(lc), height(rc))+1)
  balance(Node(avldata,lc,rc))
end

update(n::Branch{AVLData{T}}) where T = balance(n)

weigh(n::AVLNode{T}) where T = height(left(n))-height(right(n))
weigh(n::Nothing) = 0

key(d::AVLData{T}) where T = key(d.kv)
value(d::AVLData{T}) where T = value(d.kv)
 
keyval(n::AVLNode{T}) where T = data(n).kv

height(n::AVLNode{T}) where T = data(n).rank
height(n::EmptyNode) = -1

height(l::Branch{T}, r::Branch{T}) where T = max(height(l),height(r)) + 1

#!me why do I need to wrap? is that a key or a value?  should delete work on key(n)?  
#delete(n::AVLNode{T}, k::T) where T = delete(n, avldata(k))

import Base.show
show(io::IO, n::AVLData{T}) where T = print(string(n.kv,"↑",n.rank ))

