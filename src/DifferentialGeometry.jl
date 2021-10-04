module DifferentialGeometry

export Δ, poisson

using SparseArrays, LinearAlgebra

"""the cotangent of the angle between two vectors"""
cotan( v1::T, v2::T ) where {T<:AbstractArray} = (v1⋅v2)/norm(v1×v2)

"""
    Δ(F, P)

Laplace-Beltrami operator built from list of triangle indices and points.

"""
function Δ( F::Vector{T}, P::Vector{PT} ) where {T<:Union{AbstractArray,NTuple{3}}, FT, PT<:AbstractVector{FT}}
  n = length(F)*9
  indexI = Vector{Int}(undef,n)
  indexJ = Vector{Int}(undef,n)
  element = Vector{FT}(undef,n)
  s = 1
  @inbounds for ijk in F
    i,j,k = ijk
    Pᵢ = P[i]
    Pⱼ = P[j]
    Pk = P[k]
    ij = Pⱼ-Pᵢ
    jk = Pk-Pⱼ
    ki = Pᵢ-Pk
    a = cotan(ij,-ki)*FT(0.5)
    b = cotan(jk,-ij)*FT(0.5)
    c = cotan(ki,-jk)*FT(0.5)
    #L[i,i] += b+c
    indexI[s] = i; indexJ[s] = i; element[s] = b+c; s+=1
    #L[i,j] -= c
    indexI[s] = i; indexJ[s] = j; element[s] = -c; s+=1
    #L[i,k] -= b
    indexI[s] = i; indexJ[s] = k; element[s] = -b; s+=1
    #L[j,i] -= c
    indexI[s] = j; indexJ[s] = i; element[s] = -c; s+=1
    #L[j,j] += c+a
    indexI[s] = j; indexJ[s] = j; element[s] = c+a; s+=1
    #L[j,k] -= a
    indexI[s] = j; indexJ[s] = k; element[s] = -a; s+=1
    #L[k,i] -= b
    indexI[s] = k; indexJ[s] = i; element[s] = -b; s+=1
    #L[k,j] -= a
    indexI[s] = k; indexJ[s] = j; element[s] = -a; s+=1
    #L[k,k] += a+b
    indexI[s] = k; indexJ[s] = k; element[s] = a+b; s+=1
  end
  sparse(indexI, indexJ, element)
end

Δ( F::Vector{IT}, P ) where {IT<:Integer} = Δ(Iterators.partition( F, 3 ) |> collect, P)

"""
    constrain( L, pin )

Add Dirichlet boundary conditions to Linear Operator L. 
return tuple with reduced L and columns of L representing constrained dofs 
"""
function constrain_system( L_in::S, pin ) where {S<:AbstractArray}
  L = copy(L_in)
  I_ = sparse(one(eltype(S))*I,size(L)...)
  L_constrained = L[:,pin]
  L[:,pin] = I_[:,pin]
  L[pin,:] = I_[pin,:]

  (L,L_constrained)
end


"""
  add boundary conditions to rhs of linear system
"""
function constrain_rhs(rhs, u, pin, L_constrained)
  b = rhs - L_constrained*u[pin,:]

  for ic in pin
    b[ic,:] = u[ic,:]
  end
  b
end

"""
    poisson(L, b, u_boundary, boundary_index)

Solve the Poisson problem Lu = b, where L is a Laplace-beltrami operator with Dirichlet boundary specified by u_boundary at locations boundary_index.
"""
function poisson(L, b, u_boundary, boundary_index)
  (Lₒ, Lpin) = constrain_system(L, boundary_index)

  rhs = constrain_rhs(b, u_boundary, boundary_index, Lpin)

  Lₒ\rhs
end

poisson(L, b::Vector{P}, u_boundary::Vector{P}, boundary_index) where P <: AbstractVector =
  map(P, 
      eachrow(poisson(L, 
                      reduce(vcat, adjoint.(b)), 
                      reduce(vcat, adjoint.(u_boundary)), 
                      boundary_index)))

end
