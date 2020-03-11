using Test

using HalfEdges.BoundingVolumeTrees
using Base.Iterators
using AbstractTrees
using StaticArrays

const Vector3{T} = SVector{3,T}

@testset "traverse" begin
  t = Node(1, Node(-2, Node(-3),Node(3,Node(-4),nothing)), Node(2,nothing,Node(4)))
  @test length(t) == 7
  @test [data(i) for i in t] == [1,-2,-3,3,-4,2,4]
  @test [i for i in t] == reduce(vcat,DepthFirst(t))
  @test drop(DepthFirst(t),2) |> first |> data == -3
  @test drop(BreadthFirst(t),2) |> first |> data == 2
  @test data.(reduce(vcat,BreadthFirst(t))) == [1, -2, 2, -3, 3, 4, -4]
  @test search(iseven∘data,t) |> data == -2

  animaltree = reduce(insert,["dog","dawg","kawt","awnt"]; init = tree("cat"))
  @test map(data, leaves(animaltree)) == ["awnt","dawg","kawt"]
  @test tree([key(n) for n in animaltree]) == animaltree

  t = tree(["cat", "awnt", "aardvark", "ant", "dog", "dag", "zat", "doog"])
  @test leaves(t, x->length(key(x)) < 4) |> first |> key == "dag"
end

@testset "delete" begin
  t = tree(["cat", "dog","dawg","kawt","awnt", "doog", "zat", "aardvark"])
  @test map(key, leaves(t)) == ["aardvark", "dawg", "doog", "zat"]
  @test map(key, leaves(delete(t,"awnt"))) == ["aardvark", "dawg", "doog", "zat"]
  @test map(key, leaves(delete(t,"aardvark"))) == ["awnt", "dawg", "doog", "zat"]
  @test map(key, leaves(delete(t,"dog"))) == ["aardvark", "dawg", "zat"]
  @test map(key, leaves(delete(t,"cat"))) == ["aardvark", "doog", "zat"]
  @test map(key, leaves(delete(t,"diamond"))) == map(key, leaves(t))

  animals = ["cat", "ant", "emu", "dog", "zat", "aawk", "bat"]
  @test length(delete(delete(tree(animals),"cat"), "dog")) == length(tree(animals))-2
end

@testset "balanced" begin
  animals = ["cat", "dog","dawg","kawt","awnt", "doog", "zat", "aardvark"]
  @test tree((x->AVLData{String}(x,0)).(animals)) |> height == 3

  linear255 = tree(1:255|>collect)
  balanced255 = avltree(1:255|>collect)
  @test height(linear255) == 255
  @test height(balanced255) == 7

  b15 = avltree(1:15|>collect)
  sb15 = reduce((t,i)->delete(t,i),[6];init=b15)|>left|>right
  @test sb15|>height > sb15|>left|>height  
  
  bd = reduce((t,i)->delete(t,i),[12,6,2,13,4];init=b15)
  @test height(bd) == 3

  b31 = avltree(1:31|>collect)
  #seq = [16, 17, 15, 18, 5, 29, 20, 22, 17, 8, 29, 26, 11, 4]
  #@test length(bd) >= (31-14)
  #@test height(bd) == 4


  for i in 1:2000
    bd = reduce((acc,i)->delete(acc,rand(1:31)), 1:14; init = b31)
    @test length(bd) >= (31-14)
    @test height(bd) == 4
  end

  @test height(reduce((acc,_)->delete(acc,rand(1:255)), 1:127; init = balanced255)) == 7
  @test height(reduce((acc,i)->delete(acc,i), 1:129; init = balanced255)) == 6
end

@testset "traverse" begin
  t = Node(1, Node(-2, Node(-3),Node(3,Node(-4),nothing)), Node(2,nothing,Node(4)))
  @test length(t) == 7
  @test [data(i) for i in t] == [1,-2,-3,3,-4,2,4]
  @test [i for i in t] == reduce(vcat,DepthFirst(t))
  @test drop(DepthFirst(t),2) |> first |> data == -3 
  @test drop(BreadthFirst(t),2) |> first |> data == 2 
  @test data.(reduce(vcat,BreadthFirst(t))) == [1, -2, 2, -3, 3, 4, -4]
  @test search(iseven∘data,t) |> data == -2

  animaltree = reduce(insert,["dog","dawg","kawt","awnt"]; init = tree("cat")) 
  @test map(data, Leaves(animaltree)) == ["awnt","dawg","kawt"]
end


@testset "AABB" begin
  a = Vector3(1.0,1,1)
  o = Vector3(0.0,0,0)
  @test overlaps(AABB(o,a),AABB(0.5*a,2.0*a)) == true
  @test overlaps(AABB(o,a),AABB(10.0*a,20.0*a)) == false
  @test volume(AABB(o,a)) == 1.0
  @test volume(inflate(AABB(o,a),1.0)) == 3.0^3 
  @test contains(AABB(o,a),AABB(a*0.5,a*0.75)) == true
  @test contains(AABB(o,a),AABB(a,2.0*a)) == false
  @test contains(AABB(o,a),AABB(o,a)) == true
end

function aabbtree(extents, treebuilder) 
  extents = map(x->(Vector3{Float64}(x[1]),Vector3{Float64}(x[2])), extents)
  treebuilder(map( ((n,x),i)->AABBNodeData(AABB(n,x),i), extents, 1:length(extents))) 
end

@testset "AABBTree" begin
  boxes = [([-4,-1,-1],[-1,1,1]),([1,-1,-1],[4,1,1]),
           ([2,-1,-1],[3,1,1]), ([-10,-1,-1],[10,1,1]),
           ([11,-1,-1],[12,1,1]),([7,-1,-1],[8,1,1]),
           ([7.5,-1,-1],[7.75,1,1])]
  abba = aabbtree(boxes, tree)
  babba = aabbtree(boxes, avltree)
  @test contains((abba |> data |> x->x.aabb), (abba |> left |> data |> x->x.aabb))
  @test contains((abba |> data |> x->x.aabb), (abba |> right |> data |> x->x.aabb))
  @test contains((babba |> key ), (abba |> left |> key ))
  @test contains((babba |> key ), (abba |> right |> key ))
  hits = query(abba, AABB(Vector3(9.0,-10.0,-10.0),Vector3(20.0,10.0,10.0)))
  @test setdiff(hits, [4,5]) |> isempty

  t = tree(AABBNodeData(randAABB(),1))
  t = insert(t,AABBNodeData(randAABB(),2))
  t = insert(t,AABBNodeData(randAABB(),3))
  t = insert(t,AABBNodeData(randAABB(),4))
  t = insert(t,AABBNodeData(randAABB(),5))
  @test contains((t |> data |> x->x.aabb), (t |> left |> data |> x->x.aabb))
  @test contains((t |> data |> x->x.aabb), (t |> left |> right |> data |> x->x.aabb))
  @test contains((t |> left |> data |> x->x.aabb), (t |> left |> right |> data |> x->x.aabb))

end
