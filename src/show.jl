## real numbers should pritn 1.0, not suppress?

## Poly{T} is basically T[x]
## T[x] may not have an order

## to handle this case we create some functions
## which can be modified by users for other Ts

"`hasneg(::T)` attribute is true if: `pj < zero(T)` is defined."
hasneg{T}(::Type{T}) = true

"Could value possibly be negative and if so, is it?"
isneg{T}(pj::T) = hasneg(T) && pj < zero(T)

"Make `pj` positive if it is negative. (Don't call `abs` as that may not be defined, or appropriate.)"
aspos{T}(pj::T) = (hasneg(T) && isneg(pj)) ? -pj : pj

"Should a value of `one(T)` be shown as a coefficient of monomial `x^i`, `i >= 1`? (`1.0x^2` is shown, `1 x^2` is not)"
showone{T}(::Type{T}) = true


#####

## Numbers
showone{T<:Integer}(::Type{T}) = false
showone{T}(::Type{Rational{T}}) = false


## Polynomials as coefficients
hasneg{S}(::Type{Poly{S}}) = false
showoone{S}(::Type{Poly{S}}) = false



## Complex coefficients
## we say neg if real(z) < 0 || real(z) == 0 and imag(g) < 0
hasneg{T}(::Type{Complex{T}}) = true

function isneg{T}(pj::Complex{T})
    real(pj) < 0 && return true
    (real(pj) == 0 && imag(pj) < 0) && return(true)
    return false
end

showone{T}(pj::Type{Complex{T}}) = showone(T)



"Show different operations depending on mimetype. `l-` is leading minus sign."
function showop(::MIME"text/html", op)
    d = Dict("*" => "&times;", "+" => " &#43; ", "-" => " &#45; ", "l-" => "&#45;")
    d[op]
end

function showop(::MIME"text/latex", op)
    d = Dict("*" => "\\cdot ", "+" => " + ", "-" => " - ", "l-" => "-")
    d[op]
end

function showop(::MIME"text/plain", op)
    d = Dict("*" => "⋅", "+" => " + ", "-" => " - ", "l-" => "-")
    d[op]
end



###

function printpoly{T}(io::IO, p::Poly{T}, mimetype)
    first = true
    printed_anything = false
    for i in eachindex(p)
        printed = showterm(io,p,i,first, mimetype)
        first &= !printed
        printed_anything |= printed
    end
    printed_anything || print(io, zero(T))
end

function showterm{T}(io::IO,p::Poly{T},j,first, mimetype)
    pj = p[j]
    pj == zero(T) && return false

    pj = printsign(io, pj, j, first, mimetype)
    printcoefficient(io, pj, j, mimetype)
    printproductsign(io, pj, j, mimetype)
    printexponent(io,p.var,j, mimetype)
    true
end



## print the sign
## returns aspos(pj) 
function printsign{T}(io::IO, pj::T, j, first, mimetype)
    neg = isneg(pj)
    if first
        neg && print(io, showop(mimetype, "l-"))    #Prepend - if first and negative
    else
        neg ? print(io, showop(mimetype, "-")) : print(io,showop(mimetype, "+"))
    end

    aspos(pj)
end

## print * or cdot, ...
function printproductsign{T}(io::IO, pj::T, j, mimetype)
    j == 0 && return
    (showone(T) || pj != one(T)) &&  print(io, showop(mimetype, "*"))
end
  
function printcoefficient{T}(io::IO, pj::Complex{T}, j, mimetype)
        
    hasreal = abs(real(pj)) > 0
    hasimag = abs(imag(pj)) > 0
    
    if hasreal & hasimag
        print(io, '(')
        showio(io, mimetype, pj)
        print(io, ')')
    elseif hasreal
        a = real(pj)
        (showone(T) || a != one(T)) && show(io, mimetype, a)
    elseif hasimag
        b = imag(pj)
        (showone(T) || b != one(T)) && show(io,  mimetype, b)
        show(io, mimetype, im)
    else
        return
    end
end

  
## show a single term 
function printcoefficient{T}(io::IO, pj::T, j, mimetype)
    pj == one(T) && !(showone(T) || j == 0) && return
    show(io, mimetype, pj)
end

function printexponent(io,var,i, mimetype::MIME"text/latex")
    if i == 0
        return
    elseif i == 1
        print(io,var)
    else
        print(io,var,"^{$i}")
    end
end

function printexponent(io,var,i, mimetype)
    if i == 0
        return
    elseif i == 1
        print(io,var)
    else
        print(io,var,"^",i)
    end
end


####

@compat Base.show{T}(io::IO, p::Poly{T}) = show(io, MIME("text/plain"), p)
@compat function Base.show{T}(io::IO, mimetype::MIME"text/plain", p::Poly{T})
    print(io,"Poly(")
    printpoly(io, p, mimetype)
    print(io,")")

end

@compat function Base.show{T}(io::IO, mimetype::MIME"text/latex", p::Poly{T})
    print(io, "\$")
    printpoly(io, p, mimetype)
    print(io, "\$")
end

@compat function Base.show{T}(io::IO, mimetype::MIME"text/latex", a::Rational{T})
    abs(a.den) == one(T) ? print(io, a.num) : print(io, "\\frac{$(a.num)}{$(a.den)}")
end

@compat function Base.show{T<:Number}(io::IO, mimetype::MIME"text/latex", a::T)
    print(io, a)
end
