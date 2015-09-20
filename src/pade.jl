immutable Pade{T<:Number,S<:Number}
    p::Poly{T}
    q::Poly{S}
    var::Symbol
    function Pade(p::Poly{T},q::Poly{S})
        if p.var != q.var
            error("Polynomials must have same variable")
        end
        new(p, q, p.var)
    end
end
Pade{T<:Number,S<:Number}(p::Poly{T}, q::Poly{S}) = Pade{T,S}(p,q)

function Pade{T}(c::Poly{T},m::Int,n::Int)
    @assert m+n < length(c)
    rold,rnew = Poly([zeros(T,m+n+1);one(T)],c.var),Poly(c.a[1:m+n+1],c.var)
    uold,vold = Poly([one(T)],c.var),Poly([zero(T)],c.var)
    unew,vnew = vold,uold
    for i=1:n
        temp0,temp1,temp2 = rnew,unew,vnew
        q,rnew = divrem(rold,rnew)
        unew,vnew = uold-q*unew,vold-q*vnew
        rold,uold,vold = temp0,temp1,temp2

    end
    if vnew[0] == 0
        d = gcd(rnew,vnew)
        rnew /= d
        vnew /= d
    end
    Pade(rnew/vnew[0],vnew/vnew[0])
end
padeval{T}(PQ::Pade{T},x) = polyval(PQ.p,x)./polyval(PQ.q,x)
