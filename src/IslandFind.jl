"""
Fast union-find algorithm for finding disjoint islands of entities.
Uses path compression. union! and find! are not threadsafe.
"""
module IslandFind

export Archipelago, union!, find!, islands, find_islands, group_sets

struct Node{T}
  parent::Int
  rank::Int
  data::T
end

struct Archipelago
  parent::Vector{Int}
  rank::Vector{Int}
  Archipelago(n) = new(collect(1:n), fill(1, n))
end


maxentity(isls::Archipelago) = length(isl.parent)

function resize(isls::Archipelago, n)
  currn = length(isls.parent)
  for i in (currn+1):n
    push!(isls.parent, i)
    push!(isls.rank, 1)
  end
end

"""
    union!(isls, (id_a,id_b))

connect entities a and b in disjoint set data structure.
"""
function union!(isls::Archipelago, (a,b)::Tuple{Int, Int} )
  nab = max(a,b)
  if nab > length(isls.parent)
    resize(isls, nab)
  end

  roota = find!(isls, a)
  rootb = find!(isls, b)

  if roota != rootb
    ranka = isls.rank[roota]
    rankb = isls.rank[rootb]
    if ranka == rankb
      # doesn't matter who becomes parent
      isls.rank[roota] += 1
      isls.parent[rootb] = roota
    elseif ranka < rankb
      # add to higher rank parent
      isls.parent[roota] = rootb
    else
      isls.parent[rootb] = roota
    end
  end
  roota
end

"""
    find!(isls, id_a)

the representative id for the disjoint set that id_a belogs to
"""
function find!(isls::Archipelago, id::Int)
  node = id
  while true
    root = isls.parent[node]
    root == node && break
    node = root
  end
  isls.parent[id] = node
end

"""
    find_islands(connected, group_by_entity=true)

find unique disjoint set ids for each entity.  

returns different data if group_by_entity is true or false.
returns list of island ids for each entity if group_by_entity is set false
returns disjoint sets of connected entities if group_by_entity is set true

if there are entites not represented in the connections but with ids < maximum id in the connections, they will be in singleton islands.  

example:

  julia> find_islands([(1,3),(3,4)], true)
  2-element Array{Array{Int64,1},1}:
   [1, 3, 4]
   [2]

   julia> find_islands([(1,3),(3,4)], false)
   4-element Array{Int64,1}:
    1
    2
    1
    1

"""
function find_islands(connected::C, 
                      group_by_entity=true) where {E<:Union{AbstractVector, Tuple, Set}, 
                                                   C<:AbstractVector{E}}
  arp = Archipelago(first(first(connected)))
  for c in connected
    for edge in zip(c[1:end-1], c[2:end])
      union!(arp, edge)
   end
  end

  archi = map(i->find!(arp, i), 1:length(arp.parent))
  if group_by_entity
    group_sets(archi)
  else
    archi
  end
end

function find_islands(arp::Archipelago, group_by_entity=true)
  archi = map(i->find!(arp, i), 1:length(arp.parent))
  if group_by_entity
    group_sets(archi)
  else
    archi
  end
end

""" group index sets together where values at indices are equal """
function group_sets(island_id::V) where {V<:Vector{Int}}
  isls = map(i->Vector{Int}(), 1:length(island_id))  

  for (i, isl) in enumerate(island_id)
    push!(isls[isl], i)
  end
  filter(!isempty, isls)
end

end
