# Poly type manipulations

module Polynomials
#todo: sparse polynomials?

using Compat

export Poly, polyval, polyint, polyder, poly, roots
export Pade, padeval, degree

import Base: length, endof, getindex, setindex!, copy, zero, one, convert
import Base: show, print, *, /, //, -, +, ==, divrem, rem, eltype
import Base: promote_rule
if VERSION >= v"0.4"
    import Base.call
end

eps{T}(::Type{T}) = zero(T)
eps{F<:AbstractFloat}(x::Type{F}) = Base.eps(F)
eps{T}(x::Type{Complex{T}}) = eps(T)

immutable Poly{T<:Number}
    a::Vector{T}
    var::Symbol
    function Poly(a::Vector{T}, var)
        # if a == [] we replace it with a = [0]
        if length(a) == 0
            return new(zeros(T,1), symbol(var))
        else
            # determine the last nonzero element and truncate a accordingly
            a_last = max(1,findlast( p->(abs(p) > 2*eps(T)), a))
            new(a[1:a_last], symbol(var))
        end
    end
end

@compat Poly{T<:Number}(a::Vector{T}, var::Union{AbstractString,Symbol,Char}=:x) = Poly{T}(a, var)

convert{T}(::Type{Poly{T}}, p::Poly) = Poly(convert(Vector{T}, p.a), p.var)
convert{T, S<:Number}(::Type{Poly{T}}, x::S) = Poly(promote_type(T, S)[x])
convert{T, S<:Number,n}(::Type{Poly{T}}, x::Array{S,n}) = map(el->convert(Poly{promote_type(T,S)},el),x)
promote_rule{T, S}(::Type{Poly{T}}, ::Type{Poly{S}}) = Poly{promote_type(T, S)}
eltype{T}(::Poly{T}) = T

length(p::Poly) = length(p.a)
degree(p::Poly) = length(p)-1
endof(p::Poly) = length(p) - 1

getindex{T}(p::Poly{T}, i) = (i+1 > length(p.a) ? zero(T) : p.a[i+1])
function setindex!(p::Poly, v, i)
    n = length(p.a)
    if n < i+1
        resize!(p.a,i+1)
        p.a[n:i] = 0
    end
    p.a[i+1] = v
    v
end

copy(p::Poly) = Poly(copy(p.a), p.var)

zero{T}(p::Poly{T}) = Poly([zero(T)], p.var)
zero{T}(::Type{Poly{T}}) = Poly(T[])
one{T}(p::Poly{T}) = Poly([one(T)], p.var)
one{T}(::Type{Poly{T}}) = Poly([one(T)])

function show(io::IO, p::Poly)
    print(io,"Poly(")
    print(io,p)
    print(io,")")
end

function printexponent(io,var,i)
    if i == 0
    elseif i == 1
        print(io,var)
    else
        print(io,var,"^",i)
    end
end

function printterm{T}(io::IO,p::Poly{T},j,first)
    pj = p[j]
    if pj == zero(T)
        return false
    end
    neg = pj < 0
    if first
        neg && print(io, "-")    #Prepend - if first and negative
    else
        neg ? print(io, " - ") : print(io," + ")
    end
    pj = abs(pj)
    if pj != one(T) || j == 0
        show(io,pj)
    end
    printexponent(io,p.var,j)
    true
end

function printterm{T<:Complex}(io::IO,p::Poly{T},j,first)
    pj = p[j]
    abs_repj = abs(real(pj))
    abs_impj = abs(imag(pj))
    if abs_repj < 2*eps(T) && abs_impj < 2*eps(T)
        return false
    end

    # We show a negative sign either for any complex number with negative
    # real part (and then negate the immaginary part) of for complex
    # numbers that are pure imaginary with negative imaginary part

    neg = ((abs_repj > 2*eps(T)) && real(pj) < 0) ||
            ((abs_impj > 2*eps(T)) && imag(pj) < 0)

    if first
        neg && print(io, "-")    #Prepend - if first and negative
    else
        neg ? print(io," - ") : print(io," + ")
    end

    if abs_repj > 2*eps(T)    #Real part is not 0
        if abs_impj > 2*eps(T)    #Imag part is not 0
            print(io,'(',neg ? -pj : pj,')')
        else
            print(io, neg ? -real(pj) : real(pj))
        end
    else
        if abs_impj > 2*eps(T)
            print(io,'(', imag(pj),"im)")
        end
    end
    printexponent(io,p.var,j)
    true
end

function print{T}(io::IO, p::Poly{T})
    first = true
    printed_anything = false
    n = length(p)-1
    for i = 0:n
        printed = printterm(io,p,i,first)
        first &= !printed
        printed_anything |= printed
    end
    printed_anything || print(io,zero(T))
end

*{T<:Number,S}(c::T, p::Poly{S}) = Poly(c * p.a, p.var)
*{T<:Number,S}(p::Poly{S}, c::T) = Poly(p.a * c, p.var)
/(p::Poly, c::Number) = Poly(p.a / c, p.var)
-(p::Poly) = Poly(-p.a, p.var)

-(p::Poly, c::Number) = +(p, -c)
+(c::Number, p::Poly) = +(p, c)
function +(p::Poly, c::Number)
    if length(p) < 1
        return Poly([c,], p.var)
    else
        p2 = copy(p)
        p2.a[1] += c;
        return p2;
    end
end
function -{T}(c::Number, p::Poly{T})
    if length(p) < 1
        return Poly(T[c,], p.var)
    else
        p2 = -p;
        p2.a[1] += c;
        return p2;
    end
end

function +{T,S}(p1::Poly{T}, p2::Poly{S})
    if p1.var != p2.var
        error("Polynomials must have same variable")
    end
    Poly([p1[i] + p2[i] for i = 0:max(length(p1),length(p2))], p1.var)
end
function -{T,S}(p1::Poly{T}, p2::Poly{S})
    if p1.var != p2.var
        error("Polynomials must have same variable")
    end
    Poly([p1[i] - p2[i] for i = 0:max(length(p1),length(p2))], p1.var)
end


function *{T,S}(p1::Poly{T}, p2::Poly{S})
    if p1.var != p2.var
        error("Polynomials must have same variable")
    end
    R = promote_type(T,S)
    n = length(p1)-1
    m = length(p2)-1
    a = zeros(R,m+n+1)

    for i = 0:n
        for j = 0:m
            a[i+j+1] += p1[i] * p2[j]
        end
    end
    Poly(a,p1.var)
end

function divrem{T, S}(num::Poly{T}, den::Poly{S})
    if num.var != den.var
        error("Polynomials must have same variable")
    end
    m = length(den)-1
    if m == 0 && den[0] == 0
        throw(DivideError())
    end
    R = typeof(one(T)/one(S))
    n = length(num)-1
    deg = n-m+1
    if deg <= 0
        return convert(Poly{R}, zero(num)), convert(Poly{R}, num)
    end

    aQ = zeros(R, deg)
    # aR = deepcopy(num.a)
    @show num.a
    aR = R[ num.a[i] for i = 1:n+1 ]
    for i = n:-1:m
        quot = aR[i+1] / den[m]
        aQ[i-m+1] = quot
        for j = 0:m
            elem = den[j]*quot
            aR[i-(m-j)+1] -= elem
        end
    end
    pQ = Poly(aQ, num.var)
    pR = Poly(aR, num.var)

    return pQ, pR
end
/(num::Poly, den::Poly) = divrem(num, den)[1]
rem(num::Poly, den::Poly) = divrem(num, den)[2]

function ==(p1::Poly, p2::Poly)
    if p1.var != p2.var
        return false
    else
        return p1.a == p2.a
    end
end

function polyval{T,S}(p::Poly{T}, x::S)
    R = promote_type(T,S)
    lenp = length(p)
    if lenp == 0
        return zero(R)
    else
        y = convert(R, p[end])
        for i = (endof(p)-1):-1:0
            y = p[i] + x*y
        end
        return y
    end
end

polyval(p::Poly, v::AbstractVector) = map(x->polyval(p, x), v)

if VERSION >= v"0.4"
    call(p::Poly, x) = polyval(p, x)
end

function polyint{T}(p::Poly{T}, k::Number=0)
    n = length(p)
    R = typeof(one(T)/1)
    a2 = Array(R, n+1)
    a2[1] = k
    for i = 1:n
        a2[i+1] = p[i-1] / i
    end
    return Poly(a2, p.var)
end

function polyder{T}(p::Poly{T}, order::Int=1)
    n = length(p)
    if order < 0
        error("Order of derivative must be non-negative")
    elseif n <= order
        return Poly(zeros(T,0),p.var)
    elseif order == 0
        return p
    else
        a2 = Array(T, n-order)
        for i = order:n-1
            a2[i-order+1] = p[i] * prod((i-order+1):i)
        end
        return Poly(a2, p.var)
    end
end

polyint{T}(a::Array{Poly{T},1}, k::Number  = 0) = [ polyint(p,k) for p in a ]
polyder{T}(a::Array{Poly{T},1}, order::Int = 1) = [ polyder(p,order) for p in a ]
polyint{n,T}(a::Array{Poly{T},n}, k::Number  = 0) = map(p->polyint(p,k),a)
polyder{n,T}(a::Array{Poly{T},n}, order::Int = 1) = map(p->polyder(p,order),a)

# create a Poly object from its roots
function poly{T}(r::AbstractVector{T}, var=:x)
    n = length(r)
    c = zeros(T, n+1)
    c[1] = 1
    for j = 1:n
        c[2:j+1] = c[2:j+1]-r[j]*c[1:j]
    end
    return Poly(reverse(c), var)
end
poly(A::Matrix, var=:x) = poly(eig(A)[1], var)
poly(A::Matrix, var::AbstractString) = poly(eig(A)[1], symbol(var))
poly(A::Matrix, var::Char) = poly(eig(A)[1], symbol(var))

roots{T}(p::Poly{Rational{T}}) = roots(convert(Poly{promote_type(T, Float64)}, p))


# compute the roots of a polynomial
function roots{T}(p::Poly{T})
    R = promote_type(T, Float64)
    length(p) == 0 && return zeros(R, 0)

    num_leading_zeros = 0
    while abs(p[num_leading_zeros]) <= 2*eps(T)
        if num_leading_zeros == length(p)-1
            return zeros(R, 0)
        end
        num_leading_zeros += 1
    end

    num_trailing_zeros = 0
    while abs(p[end - num_trailing_zeros]) <= 2*eps(T)
        num_trailing_zeros += 1
    end

    n = endof(p)-(num_leading_zeros + num_trailing_zeros)

    n < 1 && return zeros(R, length(p) - num_trailing_zeros - 1)

    companion = zeros(R, n, n)
    an = p[end-num_trailing_zeros]
    for i = 1:n-1
        companion[i,n] = -p[num_leading_zeros + i - 1] / an
        companion[i+1,i] = 1;
    end
    companion[end,end] = -p[end-num_trailing_zeros-1] / an
    D = eigvals(companion)
    r = zeros(eltype(D),length(p)-num_trailing_zeros-1)
    r[1:n] = D
    return r
end

function gcd{T, S}(a::Poly{T}, b::Poly{S})
    #Finds the Greatest Common Denominator of two polynomials recursively using
    #Euclid's algorithm: http://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Euclid.27s_algorithm
    if all(abs(b.a).<=2*eps(S))
        return a
    else
        s, r = divrem(a, b)
        return gcd(b, r)
    end
end

include("pade.jl")

end # module Poly
