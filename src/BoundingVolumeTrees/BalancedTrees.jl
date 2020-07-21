export
AVLData,
BalancedTree,
BalanceData,
AVLTree,
AVLData,
avltree,
avlitree,
avldata,
weigh,
height,
balance,
maxheight
#build_wrapped_node

abstract type BalanceTrait end
struct Unbalanced <: BalanceTrait end
struct AVLBalanced <: BalanceTrait end
struct WAVLBalanced <: BalanceTrait end
struct RAVLBalanced <: BalanceTrait end

abstract type BalanceData{T} <: WrappedData{T} end
struct AVLData{T} <: BalanceData{T}
  kv::T
  rank::Int64
end

wrapped_data(n::AVLData{T}) where T = n.kv

AVLData(kv::T, h::S) where {T,S<:Integer} = AVLData(kv,Int64(h))
AVLData(kv::T) where T = AVLData(kv,Int64(0))
AVLData{T}(kv::T) where T = AVLData(kv,Int64(0))

const AVLNode{T} = Union{Node{AVLData{T}}, IndexedBinaryTree{AVLData{T}}}

balance_trait(::Type) = Unbalanced()
balance_trait(::Type{EmptyNode}) = Unbalanced()
balance_trait(::Type{<:AVLNode{T}}) where T = AVLBalanced()

avltree( d::T ) where T = tree( AVLData(d) )
avltree( dc::Vector{T} ) where T = tree( AVLData.(dc) )
avldata( kv::T ) where T = AVLData(kv)
avlitree( d::T ) where T = itree( AVLData(d) )
avlitree( dc::Vector{T} ) where T = itree( AVLData.(dc) )

""" maximum height of avltree with n leaves """
maxheight(n) = 1.440*log(2,n+1.065) - 0.328

function rotate_right( c::N ) where {T, N<:AVLNode{T}}
  b = left(c)
  indata = rebuild_parent_data( Val(:either), keyval(right(b)), keyval(right(c)), keyval(c) ) 
  cdata = AVLData(indata, height(right(b),right(c)))
  c = rebuild_node(Val(:left), right(b), right(c), c, cdata)

  indata = rebuild_parent_data( Val(:either), keyval(left(b)), keyval(c), keyval(b) ) 
  rebuild_node(Val(:right), left(b), c, b, AVLData(indata, height(left(b),c)))
end

function rotate_left( c::N ) where {T, N<:AVLNode{T}}
  b = right(c)
  indata = rebuild_parent_data( Val(:either), keyval(left(c)), keyval(left(b)), keyval(c) ) 
  c = rebuild_node(Val(:right), left(c), left(b), c, AVLData(indata, height(left(c),left(b))))

  indata = rebuild_parent_data( Val(:either), keyval(c), keyval(right(b)), keyval(b) ) 
  rebuild_node(Val(:left), c, right(b), b, AVLData(indata, height(c,right(b))))
end

balance( n::T ) where T = balance( n, balance_trait(T) ) 

balance( n::T, ::Unbalanced ) where T = n

function balance( st::NODE, trait::AVLBalanced ) where {T, NODE<:AVLNode{T}} 
  scale = weigh(st)
  if scale > 1
    # left heavy
    subscale = weigh(left(st))
    if subscale > 0
      # left child heavy
      rotate_right(st)
    else
      # right child heavy
      rotate_right(rebuild_node(:left, rotate_left(left(st)), right(st), st))
    end
  elseif scale < -1
    # right heavy
    subscale = weigh(right(st))
    if subscale > 0
      # left child heavy
      rotate_left(rebuild_node(:right, left(st), rotate_right(right(st)), st))
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
                      n::NODE, d::T ) where {T, NODE<:AVLNode{T}}
  avldata = AVLData{T}(d, max(height(right(n)),height(child))+1)
  balance(NODE(avldata, child, right(n)))
end

function rebuild_node(::Val{:right}, 
                      _::Branch{AVLData{T}}, 
                      child::Branch{AVLData{T}}, 
                      n::NODE, d::T ) where {T, NODE<:AVLNode{T}}
  avldata = AVLData{T}(d, max(height(left(n)),height(child))+1)
  balance(NODE(avldata,left(n),child))
end

function rebuild_node(::Val{:both}, 
                      lc::Branch{AVLData{T}}, 
                      rc::Branch{AVLData{T}}, 
                      n::NODE, d::T ) where {T, NODE<:AVLNode{T}}
  avldata = AVLData{T}(d, max(height(lc), height(rc))+1)
  balance(NODE(avldata,lc,rc))
end




function rebuild_node(dir::Val{:left}, 
                      lc::Branch{AVLData{T}}, 
                      rc::Branch{AVLData{T}}, 
                      n::IndexedBinaryTree{AVLData{T}}, d::T ) where {T}
  avldata = AVLData{T}(d, max(height(right(n)),height(lc))+1)
  balance(rebuild_node!(dir, lc, right(n), n, avldata))
end

function rebuild_node(dir::Val{:right}, 
                      lc::Branch{AVLData{T}}, 
                      rc::Branch{AVLData{T}}, 
                      n::IndexedBinaryTree{AVLData{T}}, d::T ) where {T}
  avldata = AVLData{T}(d, max(height(left(n)),height(rc))+1)
  balance(rebuild_node!(dir, left(n), rc, n, avldata))
end

function rebuild_node(dir::Val{:both}, 
                      lc::Branch{AVLData{T}}, 
                      rc::Branch{AVLData{T}}, 
                      n::IndexedBinaryTree{AVLData{T}}, d::T ) where {T}
  avldata = AVLData{T}(d, max(height(lc), height(rc))+1)
  balance(rebuild_node!(dir, lc, rc, n, avldata))
end


function ileaves(t::IndexedBinaryTree{AVLData{T}}) where T
  L = Vector{IndexedNode{AVLData{T}}}(undef, 1<<height(t))
  resize!(L, 0);
  ileaves(t.nodes, t.node, L);
  return L
end



update(n::Branch{AVLData{T}}) where T = balance(n)

weigh(n::NODE) where {T, NODE<:AVLNode{T}} = height(left(n))-height(right(n))
weigh(n::Nothing) = 0

key(d::AVLData{T}) where T = key(d.kv)
value(d::AVLData{T}) where T = value(d.kv)
 
keyval(n::NODE) where {T, NODE<:AVLNode{T}} = data(n).kv

height(n::NODE) where {T, NODE<:AVLNode{T}} = data(n).rank
height(n::EmptyNode) = -1

height(l::Branch{T}, r::Branch{T}) where T = max(height(l),height(r)) + 1

#!me why do I need to wrap? is that a key or a value?  should delete work on key(n)?  
#delete(n::AVLNode{T}, k::T) where T = delete(n, avldata(k))

import Base.show
show(io::IO, n::AVLData{T}) where T = print(string(n.kv,"↑",n.rank ))

