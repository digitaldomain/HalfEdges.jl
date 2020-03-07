# Copyright 2019 Digital Domain 3.0
#
# Licensed under the Apache License, Version 2.0 (the "Apache License")
# with the following modification; you may not use this file except in
# compliance with the Apache License and the following modification to it:
# Section 6. Trademarks. is deleted and replaced with:
#
# 6. Trademarks. This License does not grant permission to use the trade
#    names, trademarks, service marks, or product names of the Licensor
#    and its affiliates, except as required to comply with Section 4(c) of
#    the License and to reproduce the content of the NOTICE file.
#
# You may obtain a copy of the Apache License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the Apache License with the above modification is
# distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied. See the Apache License for the specific
# language governing permissions and limitations under the Apache License.

"""
A lightweight (primitive), integer-like handle to a resource.  Useful for dispatch.
i.e. pretty much an integer, but with a unique dispatch signature.
"""

export 
Handle, 
@HandleType

abstract type Handle <: Integer end

"""
    @HandleType(n)

Create a new type named n to be used as a handle to your resource.

Example:
julia> @HandleType(Foo)
Foo

julia> foo(x::Integer) = x*x
foo (generic function with 1 method)

julia> foo(Foo(2))
4

julia> foo(x::Foo) = x-x
foo (generic function with 2 methods)

julia> foo(Foo(2))
0
"""
macro HandleType(tname)
  esc(
      quote 
        primitive type $tname <: HalfEdges.Handle sizeof(Int)*8 end
        $tname(x::T) where {T<:Integer} = reinterpret($tname, Int(x))
      end
     )
end

Base.show(io::IO, x::Handle) = print(io, Int(x))

Base.ndigits0zpb( c::C, i::Integer ) where C<:Handle = ndigits0zpb( Int(c), i )
Base.promote_rule( ::Type{C}, ::Type{I} ) where {C<:Handle,I<:Integer} = C
Base.convert( ::Type{C}, number::Number ) where {C<:Handle} = C(Int(number))

Base.Int(x::Handle) = reinterpret(Int,x)
Base.UInt(x::Handle) = reinterpret(UInt,x)
Base.Int32(x::Handle) = Int32(reinterpret(Int,x))
Base.UInt32(x::Handle) = UInt32(reinterpret(UInt,x))

Base.hash(d::Handle, x::UInt64) = hash(reinterpret(UInt64,d),x)

for op in (:+, :-, :*, :/, :mod, :div, :rem, :max, :min)
  @eval Base.$op(a::H, b::H) where {H<:Handle} = H($op(Int(a),Int(b)))
end

for op in (:<, :>, :<=, :>=)
  @eval Base.$op(a::H, b::H) where {H<:Handle} = $op(Int(a),Int(b))
end

for op in (:<<, :>>)
  @eval Base.$op(a::H, b::I) where {H<:Handle, I<:Int} = H($op(Int(a),Int(b)))
end

Base.sub_with_overflow( a::H, b::H ) where {H<:Handle} = 
  ((r,f)->(H(r),f))(Base.sub_with_overflow(Int(a), Int(b))...)
Base.add_with_overflow( a::H, b::H ) where {H<:Handle} = 
  ((r,f)->(H(r),f))(Base.add_with_overflow(Int(a), Int(b))...)

