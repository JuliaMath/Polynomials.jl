# Poly type manipulations

module Polynomials
#todo: sparse polynomials?

export Poly, polyval, polyint, polyder, poly, roots

import Base: length, endof, getindex, setindex!, copy, zero, one, convert
import Base: show, print, *, /, //, -, +, ==, divrem, rem, eltype

eps{T}(::Type{T}) = zero(T)
eps{F<:FloatingPoint}(x::Type{F}) = Base.eps(F)
eps{T}(x::Type{Complex{T}}) = eps(T)

immutable Poly{T<:Number}
    a::Vector{T}
    var::Symbol
    function Poly(a::Vector{T}, var)
        new(a, symbol(var))
    end
end

Poly{T<:Number}(a::Vector{T}, var::Union(String,Symbol,Char)=:x) = Poly{T}(a, var)

convert{T}(::Type{Poly{T}}, p::Poly) = Poly(convert(Vector{T}, p.a), p.var)
promote_rule{T, S}(::Type{Poly{T}}, ::Type{Poly{S}}) = Poly{promote_type(T, S)}
eltype{T}(::Poly{T}) = T

length(p::Poly) = length(p.a)
endof(p::Poly) = length(p)

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

copy(p::Poly) = Poly(copy(p.a))

zero{T}(p::Poly{T}) = zero(Poly{T})
zero{T}(::Type{Poly{T}}) = Poly(T[])
one{T}(p::Poly{T}) = one(Poly{T})
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
    neg = abs(pj) < 0 
    if first
        neg && print(io, "-")    #Prepend - if first and negative
    else
        neg ? print(io, " - ") : print(io," + ")
    end
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
            print(io, real(pj))
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

*{T<:Number,S}(c::T, p::Poly{S}) = Poly(c * p.a)
*{T<:Number,S}(p::Poly{S}, c::T) = Poly(p.a * c)
/(p::Poly, c::Number) = Poly(p.a / c)
-(p::Poly) = Poly(-p.a)

-(p::Poly, c::Number) = +(p, -c)
+(c::Number, p::Poly) = +(p, c)
function +(p::Poly, c::Number)
    if length(p) < 1
        return Poly([c,])
    else
        p2 = copy(p)
        p2.a[1] += c;
        return p2;
    end
end
function -{T}(c::Number, p::Poly{T})
    if length(p) < 1
        return Poly(T[c,])
    else
        p2 = -p;
        p2.a[1] += c;
        return p2;
    end
end

+{T,S}(p1::Poly{T}, p2::Poly{S}) = Poly([p1[i] + p2[i] for i = 0:max(length(p1),length(p2))])
-{T,S}(p1::Poly{T}, p2::Poly{S}) = Poly([p1[i] - p2[i] for i = 0:max(length(p1),length(p2))])

function *{T,S}(p1::Poly{T}, p2::Poly{S})
    R = promote_type(T,S)
    n = length(p1)
    m = length(p2)
    a = Poly(zeros(R,m+n+2))
    for i = 0:length(p1)
        for j = 0:length(p2)
            a[i+j] += p1[i] * p2[j]
        end
    end
    a
end

function divrem{T, S}(num::Poly{T}, den::Poly{S})
    if num.var != den.var
        error("Polynomials must have same variable")
    end
    m = length(den)
    if m == 0
        throw(DivideError())
    end
    R = typeof(one(T)/one(S))
    n = length(num)
    deg = n-m+1
    if deg <= 0
        return zero(Poly{R}), convert(Poly{R}, num)
    end
    d = zeros(R, n)
    q = zeros(R, deg)
    r = zeros(R, n)
    r[:] = num.a
    for i = 1:deg
        quot = r[i] / den[1]
        q[i] = quot
        if i > 1
            d[i-1] = 0
            r[i-1] = 0
        end
        for j = 1:m
            k = i+j-1
            elem = den[j]*quot
            d[k] = elem
            r[k] -= elem
        end
    end
    return Poly(q, num.var), Poly(r, num.var)
end
/(num::Poly, den::Poly) = divrem(num, den)[1]
rem(num::Poly, den::Poly) = divrem(num, den)[2]

function ==(p1::Poly, p2::Poly)
    if p1.var != p2.var
        return false
    else
        for i = 1:max(length(p1),length(p2))
            if p1[i] != p2[i]
                return false    
            end
        end
        return true
    end
end

function polyval{T}(p::Poly{T}, x::Number)
    R = promote_type(T, typeof(x))
    lenp = length(p)
    if lenp == 0
        return zero(R)
    else
        y = convert(R, p[end])
        for i = (lenp-1):-1:0
            y = p[i] + x.*y
        end
        return y
    end
end

polyval(p::Poly, v::AbstractVector) = map(x->polyval(p, x), v)

function polyint{T}(p::Poly{T}, k::Number=0)
    n = length(p)
    R = typeof(one(T)/1)
    a2 = Array(R, n+1)
    a2[1] = k
    p2 = Poly(a2, p.var)
    for i = 1:n
        p2[i] = p[i-1] / i
    end
    p2
end

function polyder{T}(p::Poly{T})
    n = length(p)
    if n > 0
        p2 = Poly(Array(T, n-1), p.var)
        for i = 1:(n-1)
            p2[i-1] = p[i] * i
        end
        return p2
    else
        return Poly(zeros(T, 0), p.var)
    end
end

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
poly(A::Matrix, var::String) = poly(eig(A)[1], symbol(var))
poly(A::Matrix, var::Char) = poly(eig(A)[1], symbol(var))

roots{T}(p::Poly{Rational{T}}) = roots(convert(Poly{promote_type(T, Float64)}, p))

# compute the roots of a polynomial
function roots{T}(p::Poly{T})
    R = promote_type(T, Float64)
    num_zeros = 0
    if length(p) == 0
        return zeros(R, 0)
    end
    while abs(p[num_zeros]) <= 2*eps(T)
        if num_zeros == length(p)-1
            return zeros(R, 0)
        end
        num_zeros += 1
    end
    n = length(p)-num_zeros-1
    if n < 1
        return zeros(R, length(p)-1)
    end
    companion = zeros(R, n, n)
    a0 = p[num_zeros]
    for i = 1:n-1
        companion[1,i] = -p[num_zeros+i] / a0
        companion[i+1,i] = 1;
    end
    companion[1,end] = -p[end] / a0
    D,V = eig(companion)
    r = zeros(eltype(D),length(p)-1)
    r[1:n] = 1./D
    return r
end

end # module Poly