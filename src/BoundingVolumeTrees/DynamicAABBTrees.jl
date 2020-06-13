# older but much faster BVH

struct Triangles{IntT,FloatT}
  indices::Vector{IntT}
  P::Vector{SVector{3,FloatT}}
end

export 
insert, query, makeroot, isleaf, visit, visitbf, 
updateBVH, createBVH, @createBVH_nplex,
selfintersect, selfintersects, @selfintersect, @selfintersects,
DAABBNode, EmptyTree, triangle_intersect,
remove, getparent, getsibling, getleft, getright, isleft, isright, verify, getAABB

const Nullable{T} = Union{Nothing,T}
const isnull = isnothing

using Base.Iterators
#using SIMD
import Base.<

#========== macros ==========#

#== this will make a block of 
quote
  baseindex_ = 3 * (j - 1)
  j1 = A[baseindex_ + 1]
  j2 = A[baseindex_ + 2]
  j3 = A[baseindex_ + 3]
end
==#
macro simplex_indices(A,dim,j)
  Expr(:block, 
       Expr(:(=), 
            Expr(:escape, Symbol("baseindex_")),
            Expr(:call, :*, eval(dim), Expr(:call, :-, Expr(:escape, j), 1))),
       map( i->Expr(:(=), 
                    Expr(:escape, Symbol(j,i)), 
                    Expr(:ref, 
                         Expr(:escape, A), 
                         Expr(:call, :+, Expr(:escape,Symbol("baseindex_")), i))
                   ),1:eval(dim))...,
       Expr(:(=),
            Expr(:escape, Symbol(j,"vertindices")),
            Expr(:tuple, map( i->Expr(:escape, Symbol(j,i)), 1:eval(dim))...)))
end

# use to reference symbols made by other macros in same calling context
macro mref(prefix,vname)
  Expr(:escape, Symbol(prefix,vname))
end

#== 
transforms
simplex_P(P,3,i)
to 
quote
  pi1 = P[i1]
  pi2 = P[i2]
  pi3 = P[i3]
end
==# 
macro simplex_P(P,dim,j)
  Expr(:block, 
       map( i->Expr(:(=), 
                    Expr(:escape, Symbol("p",j,i)), 
                    Expr(:ref, 
                         Expr(:escape, P), 
                         Expr(:escape,Symbol(j,i)))
                   ),1:eval(dim))...)
end

macro stup(prefix,dim)
  Expr(:tuple, map(i->Expr(:escape, Symbol(prefix,i)), 1:eval(dim))...)
end

macro sargs(prefix,dim)
  Expr(:..., Expr(:tuple, map(i->Expr(:escape, Symbol(prefix,i)), 1:eval(dim))...))
end

macro sargs(prefix,iter,dim)
  Expr(:..., Expr(:tuple, map(i->Expr(:escape, Symbol(prefix,iter,i)), 1:eval(dim))...))
end

macro ifenable( isenabled, jlti, body )
  if eval(isenabled)
    :(if $(esc(jlti)) 
        $(esc(body))
      end)
  else
    esc(body)
  end
end

#============  DynamicAABBTree =============#

mutable struct DAABBNode{T,F}
  aabb::AABB{F}
  data::T
  parent::Nullable{DAABBNode{T,F}}
  left::Nullable{DAABBNode{T,F}}
  right::Nullable{DAABBNode{T,F}}
end

"""
this special type indicates a tree with nothing in it.
insert will dispatch on this type or a DAABBNode
"""
struct EmptyTree
end

const ChildNode=Nullable{DAABBNode}

DAABBNode(aabb::AABB{F}, data::T, parent::DAABBNode{T,F}, childleft::DAABBNode{T,F}, childright::DAABBNode{T,F}) where {T,F} = 
DAABBNode(aabb,data,parent,childleft,childright)

DAABBNode(aabb::AABB{F}, data::T, parent::DAABBNode{T,F}) where {T,F} = 
DAABBNode(aabb,data,parent,nothing,nothing)

DAABBNode(aabb::AABB{F}, data::T) where {T,F} = 
DAABBNode(aabb,data,nothing,nothing,nothing)

function getAABB(notree::EmptyTree)
  AABB(0.,0.,0.,0.,0.,0.)
end

function getAABB(tree::DAABBNode{T,F}) where {T,F}
  tree.aabb  
end

function getparent(kid::EmptyTree)
  EmptyTree()
end

function getparent(kid::DAABBNode{T,F}) where {T,F}
  if isnull(kid.parent)
    EmptyTree()
  else 
    kid.parent
  end
end


function getsibling(parent::DAABBNode{T,F}, kid::DAABBNode{T,F}) where {T,F}
  if !isnull(parent.right) && (parent.right == kid)
    parent.left
  else
    parent.right
  end
end

getsibling(parent::EmptyTree, kid::DAABBNode{T,F}) where {T,F} = EmptyTree()

getsibling(kid::DAABBNode{T,F}) where {T,F} = getsibling(getparent(kid),kid)

getleft(parent::EmptyTree) = EmptyTree()
getright(parent::EmptyTree) = EmptyTree()

getleft(parent::DAABBNode{T,F}) where {T,F} = isnull(parent.left) ? EmptyTree() : parent.left
getright(parent::DAABBNode{T,F}) where {T,F} = isnull(parent.right) ? EmptyTree() : parent.right

isleft(kid::DAABBNode{T,F}) where {T,F} = (getleft(getparent(kid)) == getsibling(getsibling(kid)))
isright(kid::DAABBNode{T,F}) where {T,F} = (getright(getparent(kid)) == getsibling(getsibling(kid)))

"""
return an empty tree
insert will dispatch on this type or a DAABBNode
"""
function makeroot()
  EmptyTree()
end

"""
insert a new leaf with data.  returns the new root
"""
function insert( ancestor::DAABBNode{T,F}, aabb::AABB{F}, fataabb::AABB{F}, data::T ) where {T,F}

  new_ancestor_aabb = union(fataabb,ancestor.aabb)

  if isleaf(ancestor)
    new_ancestor = DAABBNode(new_ancestor_aabb, zero(data))
    new_ancestor.parent = ancestor.parent
    new_ancestor.left = ancestor 
    new_child = DAABBNode(aabb, data, new_ancestor)
    new_ancestor.right = new_child
    ancestor.parent = new_ancestor
    new_ancestor
  else
    ancestor.aabb = new_ancestor_aabb
    vol = volume(fataabb)

    #assert( !isnull(ancestor.left) ) # we keep tree in this state. no null children.
    #assert( !isnull(ancestor.right) )
    childL = ancestor.left 
    childR = ancestor.right 

    volleft = volume(union(fataabb, childL.aabb))
    volright = volume(union(fataabb, childR.aabb))

    if volleft < volright
      childL = insert(childL, aabb, fataabb, data) 
      ancestor.left = childL
    else
      childR = insert(childR, aabb, fataabb, data) 
      ancestor.right = childR
    end
    ancestor
  end

end


function insert( ancestor::DAABBNode{T,F}, aabb::AABB{F}, data::T, slop::F=0.1 ) where {T,F}
  fataabb = inflate(aabb,slop)
  insert(ancestor, aabb, fataabb, data)
end

function insert( root::EmptyTree, aabb::AABB{F}, data::T, slop::F=0.1 ) where {T,F}
  DAABBNode(aabb, data)
end

function refitup( root::EmptyTree, shapec::ST, slop::S ) where {ST,S}
  root
end

function refitup( ancestor::DAABBNode{T,F}, shapec::ST, slop::F ) where {T,F,ST} 
  if isleaf(ancestor) 
    ancestor.aabb = AABB( shapec, ancestor.data ) 
    refitup( getparent(ancestor), shapec, slop )
  else
    #assert(!isnull(ancestor.left))
    #assert(!isnull(ancestor.right))
    childL = ancestor.left 
    childR = ancestor.right 
    laabb = childL.aabb
    raabb = childR.aabb

    if isleaf(childL)
      laabb = AABB( shapec, childL.data ) 
      childL.aabb = laabb
      laabb = inflate( laabb, slop )
    end
    if isleaf(childR)
      raabb = AABB( shapec, childR.data ) 
      childR.aabb = raabb
      raabb = inflate( raabb, slop )
    end

    ancestor.aabb = union(laabb, raabb)
    refitup( getparent(ancestor), shapec, slop )
  end
end

function remove( root::EmptyTree, parent::EmptyTree, kid::DAABBNode{T,F} ) where {T,F}
  EmptyTree()
end

function remove( root::EmptyTree, parent::DAABBNode{T,F}, kid::DAABBNode{T,F} ) where {T,F}
  # pull sibling up to parents position   
  brosis = getsibling(parent,kid)
  brosis.parent = nothing
  brosis
end

function remove( grandmapa::DAABBNode{T,F}, modad::DAABBNode{T,F}, kid::DAABBNode{T,F} ) where {T,F}
  # pull sibling up to parents position, grandparent will now raise the little bastard   
  brosis = getsibling(modad,kid)
  brosis.parent = grandmapa

  if grandmapa.left == modad
    grandmapa.left = brosis 
  else
    grandmapa.right = brosis 
  end
  grandmapa
end

function remove( root::EmptyTree, kid::DAABBNode{T,F} ) where {T,F}
  EmptyTree()
end

function remove( parent::DAABBNode{T,F}, kid::DAABBNode{T,F} ) where {T,F}
  grandparent = getparent(parent)
  remove( grandparent, parent, kid )
end

function remove( kid::DAABBNode{T,F} ) where {T,F} 
  parent = getparent(kid)
  grandparent = getparent(parent)
  remove( grandparent, parent, kid  )
end

function isleaf( child::DAABBNode{T,F} ) where {T,F} 
  isnull(child.left) && isnull(child.right)
end

#==== caching flat tree ===#

"""
flatten the tree into a vector with leaf nodes first and branch nodes after
can be used to accelerate queries
"""
function flattentree(root::DAABBNode{T,F}, nleaf::T) where {T,F}
  leafcache = Vector{DAABBNode{T,F}}(undef, nleaf)
  visit(root, node->leafcache[node.data] = node)
  branchcache = Vector{DAABBNode{T,F}}()
  nodemap = Dict{DAABBNode{T,F},T}()
  for i in 1:nleaf
    nodemap[leafcache[i]] = i
  end

  id = nleaf + 1 
  cachebranch(node) = begin
    node.data = id
    push!(branchcache,node) 
    nodemap[node] = id
    id += 1
  end
                       
  visit(root, identity, cachebranch)
  vcat(leafcache,branchcache)
end

create_childcache(flattree::Vector{DAABBNode{T,F}}) where {T,F} = map(node->( isnull(node.left) ? 0 : node.left.data, isnull(node.right) ? 0 : node.right.data), flattree);

create_aabbcache(flattree::Vector{DAABBNode{T,F}}) where {T,F} = map(node->node.aabb,flattree)

struct TreeCache{IntT,FloatT}
  aabbcache::Vector{AABB{FloatT}}
  nleaf::IntT
  childcache::Vector{Tuple{IntT,IntT}}
end

"""
create a flat memory cached verion of the tree for quick queries
"""
function cachetree( root::DAABBNode{T,F}, nleaf::T ) where {T,F}
  flattree = flattentree(root, nleaf)
  TreeCache( create_aabbcache(flattree), 
             nleaf,
             create_childcache(flattree))
end

"""
return all leaves where aabb overlaps with the queried one. uses tree cache
known bug: a single node tree will miss that hit
"""
function query( rootid::T, aabb::AABB{F}, hits::Vector{T},
                      treecache::TreeCache{T,F} ) where {T,F}
                     
  nleaf = treecache.nleaf

  nschunk = 32 
  nscapacity = nschunk
  #!me should make a new type and use that... AutoVector
  nodestack = Vector{T}(undef,nscapacity)
  inodestack = 1
  nodestack[inodestack] = rootid

  while(inodestack > 0)
    ancestor = nodestack[inodestack]
    inodestack -= 1

    if overlaps(aabb,treecache.aabbcache[ancestor])

      if ancestor <= nleaf
        push!(hits,ancestor)
      else

        childL = treecache.childcache[ancestor][1]
        childR = treecache.childcache[ancestor][2]
        if childL > 0
          inodestack += 1
          nodestack[inodestack] = childL
        end

        if childR > 0
          if inodestack == nscapacity
            nscapacity *= 2
            nodestackprev = nodestack
            nodestack = Vector{T}(undef,nscapacity)
            copyto!(nodestack, 1, nodestackprev, 1, inodestack)
          end
          inodestack += 1
          nodestack[inodestack] = childR
        end
      end
    end
  end
end

"""
return all leaves where aabb overlaps with the queried one
known bug: a single node tree will miss that hit
"""
function query( root::DAABBNode{T,F}, aabb::AABB{F}, hits::Vector{T} ) where {T,F}
  nodestack = Vector{DAABBNode{T,F}}([root])

  while(!isempty(nodestack))
    ancestor = pop!(nodestack)

    if overlaps(aabb,ancestor.aabb)

      if isleaf( ancestor )
        push!(hits,ancestor.data)
      else

        childL = ancestor.left
        childR = ancestor.right
        if !isnull(childL)
          push!(nodestack, childL)
        end

        if !isnull(childR)
          push!(nodestack, childR)
        end
      end
    end
  end
end

"""
return all aabb's in the tree the overlap the query 
"""
function query( ancestor::DAABBNode{T,F}, aabb::AABB{F} ) where {T,F}
  hits = Vector{T}()
  query( ancestor, aabb, hits )
  hits
end

"""
visit nodes in depth first order
"""
function visit( ancestor::DAABBNode{T,F}, leaffcn = show, branchfcn = identity ) where {T,F}
  if isleaf(ancestor)
    leaffcn(ancestor)
  else
    branchfcn(ancestor) 
    if !isnull(ancestor.left)
      visit( ancestor.left, leaffcn, branchfcn )
    end
    if !isnull(ancestor.right)
      visit( ancestor.right, leaffcn, branchfcn )
    end
  end
end

function visitreduce( ancestor::DAABBNode{T,F}, leaffcn, leaf0, branchfcn = (a,x)->a, branch0 = 0 ) where {T,F}
  if isleaf(ancestor)
    leaf0 = leaffcn(leaf0,ancestor)
  else
    branch0 = branchfcn(branch0,ancestor) 
    if !isnull(ancestor.left)
      leaf0, branch0 = visitreduce( ancestor.left, leaffcn, leaf0, branchfcn, branch0 )
    end
    if !isnull(ancestor.right)
      leaf0, branch0 = visitreduce( ancestor.right, leaffcn, leaf0, branchfcn, branch0 )
    end
  end
  (leaf0,branch0)
end

function visitbfq( q::Vector{DAABBNode{T,F}}, 
                        leaffcn::Function, branchfcn::Function ) where {T,F}
  if !isempty(q)
    ancestor::DAABBNode{T,F} = shift!(q)
    if isleaf(ancestor)
      leaffcn(ancestor)
    else
      branchfcn(ancestor) 
      if !isnull(ancestor.left)
        push!(q,ancestor.left)
      end
      if !isnull(ancestor.right)
        push!(q,ancestor.right)
      end
    end
    visitbfq( q, leaffcn, branchfcn )
  end
end

"""
visit each node in breadfirst order
"""
function visitbf( ancestor::DAABBNode{T,F}, leaffcn = show, branchfcn = identity ) where {T,F}
  visitbfq( [ancestor], leaffcn, branchfcn ) 
end

"""
return the depth of the tree
"""
function depth( node::DAABBNode{T,F}, curd = 1, minmax = max ) where {T,F}
  if isleaf(node)
    curd
  else
    minmax(depth(node.left, curd+1, minmax), depth( node.right, curd+1, minmax))
  end
end

#==
function printtree( q::Vector{DAABBNode{T,F}}, depth = 1, x = 1 ) where {T,F}
            if isempty(q)
               print("\n")
            else
              node = shift!(q)
              if !isleaf(node)
                  dl = depth( node.left )
                  dr = depth( node.right )
                  for i in 1:dl; print("       "); end
                  print("[  ]")
                  for i in 1:dr; print("       "); end
                  push!(q,node.left)
                  push!(q,node.right)
              else
                  if
                    ==#


"""
verify the tree satisfies bounding volume relations.  i.e. children are contained by parents
"""
function verify( node::DAABBNode{T,F}, P::Vector{Point{F}}, tris::Vector{T} ) where {T,F}
  containment = true
  if isleaf(node)
    containment = contains(node.aabb, AABB(node.data, P, tris))
  else
    if !isnull(node.left)
      containment = containment && contains(node.aabb, node.left.aabb)
      if containment
        containment = containment && verify( node.left, P, tris ) 
      end
      containment = containment && node.left.parent == node 
    end
    if !isnull(node.right)
      containment = containment && contains(node.aabb, node.right.aabb)
      if containment
        containment = containment && verify( node.right, P, tris ) 
      end
      containment = containment && node.right.parent == node 
    end
  end
  containment
end

#===== mesh stuff =====#
function createBVH( P::Vector{Point{F}}, tris::Vector{Tuple{T,T,T}}, sloppercent::F=F(0.05) ) where {F,T}
  slop = reduce(max,reduce(max,P)-reduce(min,P))*sloppercent
  root = makeroot()
  for ((i1,i2,i3), i) in zip(tris, 1:length(tris))
    root = insert(root, AABB( P[i1], P[i2], P[i3] ), i, slop) 
  end
  root
end

function createBVH( P::Vector{Point{F}}, tris::Vector{T}, sloppercent::F=F(0.05) ) where {F,T}
  slop = reduce(max,reduce(max,P)-reduce(min,P))*sloppercent
  root = makeroot()
  for ((i1,i2,i3), i) in zip(partition(tris,3), 1:div(length(tris),3))
    root = insert(root, AABB( P[i1], P[i2], P[i3] ), i, slop) 
  end
  root
end

function createBVH( trimesh::Triangles{IntT,FloatT}, sloppercent::FloatT=FloatT(0.05) ) where {IntT,FloatT}
  createBVH( trimesh.P, trimesh.indices, sloppercent )
end

macro createBVH_nplex(P, simplices, dim, sloppercent)
  quote
    slop = reduce(max,reduce(max,$(esc(P)))-reduce(min,$(esc(P))))*$(esc(sloppercent))
    root = makeroot()
    for (@stup(i,$dim), i) in zip(partition($(esc(simplices)),$dim), 1:div(length($(esc(simplices))),$dim))
      @simplex_P($(esc(P)),$dim,i)
      root = insert(root, AABB(@sargs("p",i,$dim)), i, slop) 
    end
    root
  end
end


function points(tris::Triangles{IntT,FloatT},i::IntT) where {IntT,FloatT}
  ii = (i-1)*3; 
  @inbounds r = (tris.P[tris.indices[ii+1]],tris.P[tris.indices[ii+2]],tris.P[tris.indices[ii+3]]); r
end

AABB( shapec::ST, ishape::IntT ) where {IntT<:Integer,ST} = AABB( points(shapec, ishape)... )

function updateBVH(node::DAABBNode{T,F}, 
                           shapec::ST, slop::F,
                           reinsert::Vector{DAABBNode{T,F}}) where {T,F,ST}
  #!me should compare with how box2d does it.  they do a very quick leaf only check.
  if isleaf(node)
    parent = getparent(node)
    parentAABB = getAABB(parent)
    currAABB = AABB(shapec, node.data)
    node.aabb = currAABB
    if contains(parentAABB,currAABB) 
    else
      grandmapa = getparent(parent)
      remove( grandmapa, parent, node ) 
      refitup( grandmapa, shapec, slop ) 
      #insert( root, currAABB, node.data, slop )
      push!(reinsert,node)
    end
  else
    if !isnull(node.left)
      updateBVH(  node.left , shapec, slop, reinsert )
    end
    if !isnull(node.right)
      updateBVH(  node.right , shapec, slop, reinsert )
    end
  end

end

function updateBVH( node::DAABBNode{T,F}, shapec::ST, sloppercent::F = F(0.05) ) where {T,F,ST}
  #return createBVH( shapec )

  # create a foster grandparent to deal with depth 1 child removal ( will promote new root )
  knownuniverse = AABB(-Inf,Inf,-Inf,Inf,-Inf,Inf)
  fostergparent = DAABBNode( knownuniverse, 0 )
  fostergparent.left = node;

  # need well-formed children, so create a valid, but irrelevant child that fits in AABB of main node
  fostergparent.right = DAABBNode( AABB(shapec,1), 1 )
  node.parent = fostergparent

  P = shapec.P
  slop = reduce(max,reduce(max,P)-reduce(min,P))*sloppercent
  reinsert = Vector{DAABBNode{T,F}}()
  updateBVH( node, shapec, slop, reinsert ) 

  root = getleft(fostergparent)
  if (typeof(root) != EmptyTree) 
    root.parent = nothing
  end

  #println("verify after refit/remove: ", verify(root, shapec))
  for inode in reinsert
     root = insert( root, inode.aabb, inode.data, slop )
  end
  #println("verify after reinsert: ", verify(root, shapec))
  root
end

macro selfintersect( BVH, Pa, indsa, dima, Pquery, indsq, dimq, collidef, skipduplicates )
  quote
    nleaf = div(length($(esc(indsa))),$dima)
    tc = cachetree($(esc(BVH)),nleaf)
    rootid = $(esc(BVH)).data

    nquery = div(length($(esc(indsq))),$dimq)
    hits = map(x->Vector{Tuple{eltype($(esc(indsa))),eltype($(esc(indsa))),Any}}(),1:nquery)
    Threads.@threads for i in 1:nquery
      @simplex_indices($(esc(indsq)),$dimq,i)
      @simplex_P($(esc(Pquery)),$dimq,i)

      elhits = Vector{typeof(rootid)}() 

      query(rootid, AABB(@sargs("p",i,$dimq)), elhits, tc)

      for j in elhits
        # we can only skip duplicates if we are collding same simplex types on same mesh
        @ifenable $skipduplicates j < i begin
          @simplex_indices($(esc(indsa)),$dima,j)

          #!me could use @ifenable here to select between selfintersection and non self
          if isempty(intersect(@mref(i,"vertindices"),@mref(j,"vertindices")))
            @simplex_P($(esc(Pa)),$dima,j)

            #println("BVH check ", j, ", ", i)
            #hitfound = $(esc(collidef))(@sargs("p",i,$dimq),@sargs("p",j,$dima))
            #println("check hit: ", (@sargs("p",j,$dima), @sargs("p",i,$dimq))) #!me
            hit = $(esc(collidef))(@sargs("p",j,$dima), @sargs("p",i,$dimq))
            if first(hit)
              push!(hits[i],(j,i,hit)) 
            end
          end
        end
      end
    end
    reduce(vcat,hits)
  end
end

function selfintersect( BVH::DAABBNode{T,F}, triangles::Triangles{T,F}, tricollide::FN ) where {T,F,FN}
  P = triangles.P
  tris = triangles.indices
  @selfintersect( BVH, P, tris, 3, P, tris, 3, tricollide, true )
end

macro selfintersects( BVH, Pa, indsa, dima, Pquery, indsq, dimq, collidef, skipduplicates )
  quote
    nleaf = div(length($(esc(indsa))),$dima)
    tc = cachetree($(esc(BVH)),nleaf)
    rootid = $(esc(BVH)).data

    nquery = div(length($(esc(indsq))),$dimq)
    #for i in 1:ntri #nothreads
    hashit = map(x->false,1:nquery) #threads
    Threads.@threads for i in 1:nquery #threads
      @simplex_indices($(esc(indsq)),$dimq,i)
      @simplex_P($(esc(Pquery)),$dimq,i)

      elhits = Vector{typeof(rootid)}() 

      #query(BVH, AABB(p1, p2, p3), elhits)
      query(rootid, AABB(@sargs("p",i,$dimq)), elhits, tc) #threads

      for j in elhits
        #!me ack, this optimization only works if we are colliding same types... otherwise, super wrong.
        #!me can selectively add with another macro pass...
        #if j < i
        @ifenable $skipduplicates j < i begin
          @simplex_indices($(esc(indsa)),$dima,j)

          if isempty(intersect(@mref(i,"vertindices"),@mref(j,"vertindices")))
            @simplex_P($(esc(Pa)),$dima,j)

            #hitfound = triangle_intersect(pi1,pi2,pi3, pj1, pj2, pj3,cfunc) #nothreads 
            #hitfound = $(esc(collidef))(@sargs("p",i,$dimq),@sargs("p",j,$dima)) #threads 
            hitfound = $(esc(collidef))(@sargs("p",j,$dima), @sargs("p",i,$dimq)) #threads 
            #println("test:", @sargs("p",i,$dimq),@sargs("p",j,$dima), hitfound)
            if first(hitfound)
              hashit[i] = true 
              #return true #nothreads
            end
          end
        end
      end
    end
    #return false #nothreads
    return reduce((accum,el)->accum || el,hashit; init=false) #threads
  end
end

function selfintersects( BVH::DAABBNode{T,F}, triangles::Triangles{T,F}, tricollide::FN ) where {T,F,FN}
  P = triangles.P
  tris = triangles.indices
  @selfintersects( BVH, P, tris, 3, P, tris, 3, tricollide, false )
end

function selfintersectdumb( P::Vector{Point{F}}, tris::Vector{T}, trintersect::FN ) where {T,F,FN}
  hits = Vector{Tuple{Int64,Int64,Any}}()
  for (vertindices, i) in zip(partition(tris,3), 1:div(length(tris),3))
    (i1,i2,i3) = vertindices
    p1 = P[i1]
    p2 = P[i2]
    p3 = P[i3]
    paabb = AABB(p1, p2, p3)

    for (vertindicesj, j) in zip(partition(tris,3), 1:div(length(tris),3))
      if i < j
        (j1,j2,j3) = vertindicesj
        p1j = P[j1]
        p2j = P[j2]
        p3j = P[j3]
        paabbj = AABB(p1j, p2j, p3j)
        if isempty(intersect(vertindices,(j1,j2,j3))) && overlaps(paabb,paabbj)
          hitfound = trintersect(p1,p2,p3, p1j,p2j,p3j) 
          if first(hitfound)
            push!(hits,(i,j,hitfound))
          end
        end
      end
    end
  end

  unique(hits)
end
