# Poly type manipulations

module Polynomial
#todo: division
#todo: sparse polynomials?

export Poly, polyval, polyint, polydir, poly, roots

import Base.length, Base.endof, Base.ref, Base.assign, Base.copy, Base.zero, Base.one
import Base.show, Base.*, Base./, Base.-, Base.+, Base.==

type Poly{T<:Number}
    a::Vector{T}
    nzfirst::Int #for effiencicy, track the first non-zero index
    function Poly(a::Vector{T})
        nzfirst = 0 #find and chop leading zeros
        for i = 1:length(a)
            if a[i] != 0 then
                break
            end
            nzfirst = i
        end
        new(a, nzfirst)
    end
end

Poly{T<:Number}(a::Vector{T}) = Poly{T}(a)

length(p::Poly) = length(p.a)-p.nzfirst
endof(p::Poly) = length(p)
ref(p::Poly, i) = p.a[i+p.nzfirst]
assign(p::Poly, v, i) = (p.a[i+p.nzfirst] = v)

copy(p::Poly) = Poly(copy(p.a[1+p.nzfirst:end]))

zero{T}(p::Poly{T}) = Poly([zero(T)])
one{T}(p::Poly{T}) = Poly([one(T)])

function show(io::IO,p::Poly)
    n = length(p)
    print(io,"Poly(")
    if n <= 0
        print(io,"0")
    elseif n == 1
        print(io,p[1])
    else
        print(io,"$(p[1])x^$(n-1)");
        for i = 2:n-1
            if p[i] != 0
                print(io," + $(p[i])x^$(n-i)")
            end
        end
        if p[n] != 0
            print(io," + $(p[n])")
        end
    end
    print(io,")")
end

function show{T<:Complex}(p::Poly{T})
    n = length(p)
    print("Poly(")
    if n <= 0
        print("0")
    elseif n == 1
        print("[$(p[1])]")
    else
        print("[$(p[1])]x^$(n-1)")
        for i = 2:n-1
            if p[i] != 0
                print(" + [$(p[i])]x^$(n-i)")
            end
        end
        if p[n] != 0
            print(" + [$(p[n])]")
        end
    end
    print(")")
end

*(c::Number, p::Poly) = Poly(c * p.a[1+p.nzfirst:end])
*(p::Poly, c::Number) = Poly(c * p.a[1+p.nzfirst:end])
/(p::Poly, c::Number) = Poly(p.a[1+p.nzfirst:end] / c)
-(p::Poly) = Poly(-p.a[1+p.nzfirst:end])

-(p::Poly, c::Number) = +(p, -c)
+(c::Number, p::Poly) = +(p, c)
function +(p::Poly, c::Number)
    if length(p) < 1
        return Poly([c,])
    else
        p2 = copy(p)
        p2.a[end] += c;
        return p2;
    end
end
function -(c::Number, p::Poly)
    if length(p) < 1
        return Poly([c,])
    else
        p2 = -p;
        p2.a[end] += c;
        return p2;
    end
end

function +{T,S}(p1::Poly{T}, p2::Poly{S})
    R = promote_type(T,S)
    n = length(p1)
    m = length(p2)
    if n > m
        a = Array(R, n)
        for i = 1:m
            a[n-m+i] = p1[n-m+i] + p2[i]
        end
        for i = 1:n-m
            a[i] = p1[i]
        end
    else
        a = Array(R, m)
        for i = 1:n
            a[m-n+i] = p1[i] + p2[m-n+i]
        end
        for i = 1:m-n
            a[i] = p2[i]
        end
    end
    Poly(a)
end

function -{T,S}(p1::Poly{T}, p2::Poly{S})
    R = promote_type(T,S)
    n = length(p1)
    m = length(p2)
    if n > m
        a = Array(R, n)
        for i = 1:m
            a[n-m+i] = p1[n-m+i] - p2[i]
        end
        for i = 1:n-m
            a[i] = p1[i]
        end
    else
        a = Array(R, m)
        for i = 1:n
            a[m-n+i] = p1[i] - p2[m-n+i]
        end
        for i = 1:m-n
            a[i] = -p2[i]
        end
    end
    Poly(a)
end

function *{T,S}(p1::Poly{T}, p2::Poly{S})
    R = promote_type(T,S)
    n = length(p1)
    m = length(p2)
    if n == 0 || m == 0
        return Poly(R[])
    end
    a = zeros(R, n+m-1)
    for i = 1:length(p1)
        for j = 1:length(p2)
            a[i+j-1] += p1[i] * p2[j]
        end
    end
    Poly(a)
end


function ==(p1::Poly, p2::Poly)
    if length(p1) != length(p2)
        return false
    else
        return p1.a[1+p1.nzfirst:end] == p2.a[1+p2.nzfirst:end]
    end
end

function polyval{T}(p::Poly{T}, x::Number)
    R = promote_type(T, typeof(x))
    lenp = length(p)
    if lenp == 0
        return zero(R)
    else
        y = convert(R, p[1])
        for i = 2:lenp
            y = p[i] + x.*y
        end
        return y
    end
end

function polyval(p::Poly, x::AbstractVector)
    y = zeros(size(x))
    for i = 1:length(x)
        y[i] = polyval(p, x[i])
    end
    return y
end

polyint(p::Poly) = polyint(p, 0)
function polyint{T}(p::Poly{T}, k::Number)
    R = promote_type(promote_type(T, Float64), typeof(k))
    n = length(p)
    a2 = Array(R, n+1)
    for i = 1:n
        a2[i] = p[i] / (n-i+1)
    end
    a2[end] = k
    Poly(a2)
end

function polydir{T}(p::Poly{T})
    n = length(p)
    if n > 0
        a2 = Array(T, n-1)
        for i = 1:n-1
            a2[i] = p[i] * (n-i)
        end
    else
        a2 = zeros(T, 0)
    end
    Poly(a2)
end

# create a Poly object from its roots
function poly{T}(r::AbstractVector{T})
    n = length(r)
    c = zeros(T, n+1)
    c[1] = 1
    for j = 1:n
        c[2:j+1] = c[2:j+1]-r[j]*c[1:j]
    end
    return Poly(c)
end
poly(A::Matrix) = poly(eig(A)[1])

# compute the roots of a polynomial
function roots{T}(p::Poly{T})
    num_zeros = 0
    if length(p) == 0 return zeros(T,0) end
    while p[end-num_zeros] == 0
        if num_zeros == length(p)-1
            return zeros(T, 0)
        end
        num_zeros += 1
    end
    n = length(p)-num_zeros-1
    if n < 1
        return zeros(T, length(p)-1)
    end
    R = promote_type(T, Float64)
    companion = zeros(R, n, n)
    a0 = p[end-num_zeros]
    for i = 1:n-1
        companion[1,i] = -p[end-num_zeros-i] / a0
        companion[i+1,i] = 1;
    end
    companion[1,end] = -p[1] / a0
    D,V = eig(companion)
    T_r = typeof(real(D[1]))
    T_i = typeof(imag(D[1]))
    if all(imag(D) .< 2*eps(T_i))
        r = zeros(T_r, length(p)-1)
        r[1:n] = 1./real(D)
        return r
    else
        r = zeros(typeof(D[1]),length(p)-1)
        r[1:n] = 1./D
        return r
    end
end

end # module Poly
