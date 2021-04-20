## An example subtype of AbstractRationalFunction
## that might prove useful as a guide
module TransferFunction
# https://github.com/andreasvarga/DescriptorSystems.jl/blob/main/src/types/RationalFunction.jl

using Polynomials
import Polynomials: AbstractPolynomial, AbstractRationalFunction, RationalFunction
import Polynomials: pqs

export RationalTransferFunction
export sampling_time, gain

struct RationalTransferFunction{T,X,P<:AbstractPolynomial{T,X},Ts} <: AbstractRationalFunction{T,X,P}
    num::P
    den::P
    function RationalTransferFunction{T,X,P,Ts}(num::P, den::P) where{T,X,P<:AbstractPolynomial{T,X}, Ts}
        check_den(den)
        new{T,X,P,Ts}(num,den)
    end
    function RationalTransferFunction{T,X,P,Ts}(num::P, den::P,ts::Union{Real,Nothing}) where{T,X,P<:AbstractPolynomial{T,X}, Ts}
        check_den(den)        
        check_Ts(Ts,ts)        
        new{T,X,P,Ts}(num,den)
    end
    # can promote constants to polynomials too
    function  RationalTransferFunction{T,X,P,Ts}(num::S, den::P, ts::Union{Real,Nothing}) where{S, T,X,P<:AbstractPolynomial{T,X}, Ts}
        check_den(den)
        check_Ts(Ts,ts)        
        new{T,X,P,Ts}(convert(P,num),den)
    end
    function  RationalTransferFunction{T,X,P,Ts}(num::P,den::S, ts::Union{Real,Nothing}) where{S, T,X,P<:AbstractPolynomial{T,X}, Ts}
        check_den(den)
        check_Ts(Ts,ts)
        new{T,X,P,Ts}(num, convert(P,den))
    end
    function RationalTransferFunction{T,X,P}(num::P, den::P, Ts::Union{Real,Nothing}) where{T,X,P<:AbstractPolynomial{T,X}}

        check_den(den)
        Ts′ = standardize_Ts(Ts)
        new{T,X,P,Val(Ts′)}(num,den)
    end
end

(pq::RationalTransferFunction)(x) = Polynomials.eval_rationalfunction(x, pq)

function Base.convert(::Type{PQ}, pq::RationalTransferFunction) where {PQ <:RationalFunction}
    p,q = pq
    p//q
end

# alternate constructor
function RationalTransferFunction(p′::P, q′::Q, Ts::Union{Real,Nothing}) where {T,X,P<:AbstractPolynomial{T,X},
                                                                                S,  Q<:AbstractPolynomial{S,X}}

    p,q = promote(p′, q′)
    R = eltype(p)
    RationalTransferFunction{R,X,typeof(p)}(p,q,Ts)
end


function Polynomials.rational_function(::Type{PQ}, p::P, q::Q) where {PQ <:RationalTransferFunction,
                                                          T,X,   P<:AbstractPolynomial{T,X},
                                                          S,   Q<:AbstractPolynomial{S,X}}
    RationalTransferFunction(promote(p,q)..., sampling_time(PQ))
end



## helpers for constructors
# standardize Ts or throw error
function standardize_Ts(Ts)
    isnothing(Ts) || Ts >= 0 || Ts == -1 || 
        throw(ArgumentError("Ts must be either a positive number, 0 (continuous system), or -1 (unspecified)"))
    Ts′ = isnothing(Ts) ? Ts : Float64(Ts)
end
function check_Ts(Ts, ts)
    ValT(Ts) == promote_Ts(ValT(Ts), ts) || throw(ArgumentError("sampling times have mismatch"))
end
function check_den(den)
    iszero(den) && throw(ArgumentError("Cannot create a rational function with zero denominator"))
end

ValT(::Val{T}) where {T} = T
sampling_time(pq::RationalTransferFunction{T,X,P,Ts}) where {T,X,P,Ts} = ValT(Ts)
sampling_time(::Type{𝑷}) where {T,X,P,Ts, 𝑷<:RationalTransferFunction{T,X,P,Ts}} = ValT(Ts)



## ----

function Base.convert(PQ′::Type{PQ}, p::P) where {PQ <: RationalTransferFunction, P<:AbstractPolynomial}
    PQ(p, one(p), sampling_time(PQ))
end
function Base.convert(PQ′::Type{PQ}, p::Number) where {PQ <: RationalTransferFunction}
    PQ(p, one(eltype(PQ)), sampling_time(PQ))
end

function promote_Ts(p,q)
    Ts1,Ts2 = sampling_time.((p,q))
    promote_Ts(Ts1, Ts2)
end

function promote_Ts(Ts1::Union{Float64,Nothing}, Ts2::Union{Float64,Nothing})
    isnothing(Ts1) && (return Ts2)
    isnothing(Ts2) && (return Ts1)
    Ts1 == Ts2 && (return Ts1)  
    Ts1 == -1 && (Ts2 > 0 ? (return Ts2) : error("Sampling time mismatch"))
    Ts2 == -1 && (Ts1 > 0 ? (return Ts1) : error("Sampling time mismatch"))
    error("Sampling time mismatch")
end

function Base.promote_rule(::Type{PQ}, ::Type{PQ′}) where {T,X,P,Ts,PQ <: RationalTransferFunction{T,X,P,Ts},
                                                           T′,X′,P′,Ts′,PQ′ <: RationalTransferFunction{T′,X′,P′,Ts′}}

    S = promote_type(T,T′)
    Polynomials.assert_same_variable(X,X′)
    Y = X
    Q = promote_type(P, P′)
    ts = promote_Ts(PQ, PQ′)
    RationalTransferFunction{S,Y,Q,Val(ts)}
end
                    
                         



# zeros (roots) and poles are defined in common. The gain is specific to this  type
"""
    gain(pq::RationalTransferFunction)

The ratio of the leading coefficients
"""
function gain(pq::PQ) where {PQ <: RationalTransferFunction}
    p,q = pqs(pq)
    p[end]/q[end]
end

    
# need to check here
#
"""
     rt = adjoint(r)
Compute the adjoint `rt(λ)` of the rational transfer function `r(λ)` such that for 
`r(λ) = num(λ)/den(λ)` we have:
    (1) `rt(λ) = conj(num(-λ))/conj(num(-λ))`, if `r.Ts = 0`; 
    (2) `rt(λ) = conj(num(1/λ))/conj(num(1/λ))`, if `r.Ts = -1` or `r.Ts > 0`.
"""
function Base.adjoint(pq::PQ) where {PQ <: RationalTransferFunction}
    p,q = pqs(pq)
    Ts = sampling_time(pq)
    if Ts != nothing && iszero(Ts)
        # p(-λ)/q(-λ)
        p′ = poly_scale(p, -1)
        q′ = poly_scale(q, -1)
        return RationalTransferFunction(p′, q′, Ts)
    else
       # p(1/λ) / q(1/λ) = poly_inversion(p) / poly_inversion(q)
       # maps oo -> 0
       p′ = poly_inversion(p)
       q′ = poly_inversion(q)
       return RationalTransferFunction(p′, q′, Ts)
    end
end

## XXX confmap ...



end
