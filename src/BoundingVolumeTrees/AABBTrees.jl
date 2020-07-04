export
AABBNode,
AABBNodeData,
entity,
query

import Base.isless
using Base.Iterators

struct AABBNodeData{T,S}
  aabb::AABB{T}
  entity::S
end

AABBNodeData{T,S}() where {T,S} = AABBNodeData{T,S}(AABB{T}(), nothing)

const AABBNode{T,S} = Node{AABBNodeData{T,S}}

# Traits #
abstract type QueryTrait  end
struct IntervalTree <: QueryTrait end
struct SearchTree <: QueryTrait end

key(n::AABBNodeData{T,S}) where {T,S} = n.aabb 

fullness(::Type{Node{AABBNodeData{T,S}}}) where {T,S} = FullTree

query_trait(::Type) = SearchTree()
query_trait(::Type{<:Node{AABBNodeData{T,S}}}) where {T,S} = IntervalTree()
query_trait(::Type{<:Node{<:BalanceData{AABBNodeData{T,S}}}}) where {T,S} = IntervalTree()

entity(n::AABBNode{T,S}) where {T,S} = data(n).entity
entity(n::AABBNodeData{T,S}) where {T,S} = n.entity
entity(_::Type{S}, n::AABBNodeData{T,S}) where{T,S} = n.entity 
entity(s::Type{S}, n::Nothing) where {S} = zero(S) 

function branch_compare( aabb::AABB{T} ) where {T}
  n -> begin
    if left(n) == nothing || right(n) == nothing ||
      volume(union(aabb,key(left(n)))) < volume(union(aabb,key(right(n))))
      :left
    else
      :right
    end
  end
end

partial(f,x) = y->f(x,y)

function rebuild_parent_data(direction, 
                             lcd::AABBNodeData{T,S}, 
                             rcd::AABBNodeData{T,S}, 
                             parent::Union{Nothing, AABBNodeData{T,S}}) where {T,S}
  larger_aabb = union(lcd.aabb,rcd.aabb)
  AABBNodeData{T,S}(larger_aabb, entity(S,parent))
end

function rebuild_parent_data(direction, 
                             lcd::Nothing, 
                             rcd::AABBNodeData{T,S}, 
                             parent::Union{Nothing, AABBNodeData{T,S}}) where {T,S}
  AABBNodeData{T,S}(rcd.aabb, entity(S,parent))
end

function rebuild_parent_data(direction, 
                             lcd::AABBNodeData{T,S}, 
                             rcd::Nothing, 
                             parent::Union{Nothing, AABBNodeData{T,S}}) where {T,S}
  AABBNodeData{T,S}(lcd.aabb, entity(S,parent))
end

function rebuild_parent_data(direction, 
                             lcd::Nothing, 
                             rcd::Nothing, 
                             parent::AABBNodeData{T,S}) where {T,S}
  parent 
end

insert( t::AABBNode{T,S}, aabb::AABB{T}, entity::S ) where {T,S} = insert( t, AABBNodeData(aabb, entity) )

query( t::AABBNode{T,S}, qaabb::AABB{T} ) where {T,S} = entity.(Leaves(t,n->overlaps(data(n).aabb, qaabb )))

