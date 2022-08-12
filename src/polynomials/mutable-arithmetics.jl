using MutableArithmetics
const MA = MutableArithmetics

function _resize_zeros!(v::Vector, new_len)
    old_len = length(v)
    if old_len < new_len
        resize!(v, new_len)
        for i in (old_len + 1):new_len
            v[i] = zero(eltype(v))
        end
    end
end

"""
    add_conv(out::Vector{T}, E::Vector{T}, k::Vector{T})
Returns the vector `out + fastconv(E, k)`. Note that only
`MA.buffered_operate!` is implemented.
"""
function add_conv end

# The buffer we need is the buffer needed by the `MA.add_mul` operation.
# For instance, `BigInt`s need a `BigInt` buffer to store `E[x] * k[i]` before
# adding it to `out[j]`.
function MA.buffer_for(::typeof(add_conv), ::Type{Vector{T}}, ::Type{Vector{T}}, ::Type{Vector{T}}) where {T}
    return MA.buffer_for(MA.add_mul, T, T, T)
end
function MA.buffered_operate!(buffer, ::typeof(add_conv), out::Vector{T}, E::Vector{T}, k::Vector{T}) where {T}
    for x in eachindex(E)
        for i in eachindex(k)
            j = x + i - 1
            out[j] = MA.buffered_operate!(buffer, MA.add_mul, out[j], E[x], k[i])
        end
    end
    return out
end

"""
    @register_mutable_arithmetic
Register polynomial type (with vector based backend) to work with MutableArithmetics
"""
macro register_mutable_arithmetic(name)
    poly = esc(name)
    quote
        MA.mutability(::Type{<:$poly}) = MA.IsMutable()

        function MA.promote_operation(::Union{typeof(+), typeof(*)},
                                      ::Type{$poly{S,X}}, ::Type{$poly{T,X}}) where {S,T,X}
            R = promote_type(S,T)
            return $poly{R,X}
        end

        function MA.buffer_for(::typeof(MA.add_mul),
                               ::Type{<:$poly{T,X}},
                               ::Type{<:$poly{T,X}}, ::Type{<:$poly{T,X}}) where {T,X}
            V = Vector{T}
            return MA.buffer_for(add_conv, V, V, V)
        end

        function MA.buffered_operate!(buffer, ::typeof(MA.add_mul),
                                              p::$poly, q::$poly, r::$poly)
            ps, qs, rs = coeffs(p), coeffs(q), coeffs(r)
            _resize_zeros!(ps, length(qs) + length(rs) - 1)
            MA.buffered_operate!(buffer, add_conv, ps, qs, rs)
            return p
        end
    end
end

## Ambiguities. Issue #435
#Base.:+(p::P, ::MutableArithmetics.Zero) where {T, X, P<:Polynomials.AbstractPolynomial{T, X}} = p
#Base.:+(p::P, ::T) where {T<:MutableArithmetics.Zero, P<:Polynomials.StandardBasisPolynomial{T}} = p
