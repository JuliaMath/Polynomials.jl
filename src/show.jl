## Poly{T} is basically T[x], with T a Ring.
## T[x] may not have an order so abs, comparing to 0 may not be defined.

## to handle this case we create some functions
## which can be modified by users for other Ts

"`hasneg(::T)` attribute is true if: `pj < zero(T)` is defined."
hasneg(::Type{T}) where {T} = false

"Could value possibly be negative and if so, is it?"
isneg(pj::T) where {T} = hasneg(T) && pj < zero(T)

"Make `pj` positive if it is negative. (Don't call `abs` as that may not be defined, or appropriate.)"
aspos(pj::T) where {T} = (hasneg(T) && isneg(pj)) ? -pj : pj

"Should a value of `one(T)` be shown as a coefficient of monomial `x^i`, `i >= 1`? (`1.0x^2` is shown, `1 x^2` is not)"
showone(::Type{T}) where {T} = true


#####

## Numbers
hasneg(::Type{T}) where {T<:Real} = true

### Integer
showone(::Type{T}) where {T<:Integer} = false
showone(::Type{Rational{T}}) where {T<:Integer} = false




### Complex coefficients
hasneg(::Type{Complex{T}}) where {T} = true      ## we say neg if real(z) < 0 || real(z) == 0 and imag(g) < 0

function isneg(pj::Complex{T}) where {T}
    real(pj) < 0 && return true
    (real(pj) == 0 && imag(pj) < 0) && return(true)
    return false
end

showone(pj::Type{Complex{T}}) where {T} = showone(T)


### Polynomials as coefficients
hasneg(::Type{Poly{S}}) where {S} = false
showone(::Type{Poly{S}}) where {S} = false



#####

"Show different operations depending on mimetype. `l-` is leading minus sign."
function showop(::MIME"text/plain", op)
    d = Dict("*" => "*", "+" => " + ", "-" => " - ", "l-" => "-")
    d[op]
end

function showop(::MIME"text/latex", op)
    d = Dict("*" => "\\cdot ", "+" => " + ", "-" => " - ", "l-" => "-")
    d[op]
end

function showop(::MIME"text/html", op)
    d = Dict("*" => "&#8729;", "+" => " &#43; ", "-" => " &#45; ", "l-" => "&#45;")
    d[op]
end



###

"""
    printpoly(io::IO, p::Poly, mimetype = MIME"text/plain"(); descending_powers=false, offset::Int=0)

Print a human-readable representation of the polynomial `p` to `io`. The MIME
types "text/plain" (default), "text/latex", and "text/html" are supported. By
default, the terms are in order of ascending powers, matching the order in
`coeffs(p)`; specifying `descending_powers=true` reverses the order.
`offset` allows for an integer number to be added to the exponent, just for printing.

# Examples
```jldoctest
julia> printpoly(stdout, Poly([1,2,3], :y))
1 + 2*y + 3*y^2
julia> printpoly(stdout, Poly([1,2,3], :y), descending_powers=true)
3*y^2 + 2*y + 1
julia> printpoly(stdout, Poly([2, 3, 1], :z), descending_powers=true, offset=-2)
1 + 3*z^-1 + 2*z^-2
julia> printpoly(stdout, Poly([-1, 0, 1], :z), offset=-1, descending_powers=true)
z - z^-1
```
"""
function printpoly(io::IO, p::Poly{T}, mimetype=MIME"text/plain"(); descending_powers=false, offset::Int=0) where {T}
    first = true
    printed_anything = false
    for i in (descending_powers ? reverse(eachindex(p)) : eachindex(p))
        printed = showterm(io, p[i], p.var, i+offset, first, mimetype)
        first &= !printed
        printed_anything |= printed
    end
    printed_anything || print(io, zero(T))
    return nothing
end

"""
    showterm(io::IO, pj, var, j, first, mimetype)

Show the term `pj * var^j`.
Returns `true` after successfully printing.
"""
function showterm(io::IO, pj::T, var, j, first::Bool, mimetype) where {T}
    pj == zero(T) && return false

    pj = printsign(io, pj, first, mimetype)
    !(pj == one(T) && !(showone(T) || j == 0)) && printcoefficient(io, pj, j, mimetype)
    printproductsign(io, pj, j, mimetype)
    printexponent(io, var, j, mimetype)
    true
end



## print the sign
## returns aspos(pj)
function printsign(io::IO, pj::T, first, mimetype) where {T}
    neg = isneg(pj)
    if first
        neg && print(io, showop(mimetype, "l-"))    #Prepend - if first and negative
    else
        neg ? print(io, showop(mimetype, "-")) : print(io,showop(mimetype, "+"))
    end

    aspos(pj)
end

## print * or cdot, ...
function printproductsign(io::IO, pj::T, j, mimetype) where {T}
    j == 0 && return
    (showone(T) || pj != one(T)) &&  print(io, showop(mimetype, "*"))
end

# show a single term
# Other types can overload Polynomials.printcofficient with a mimetype
# or Base.show_unquoted(io, pj, indent, prec)
# For example

"""
    printcoefficient(io::IO, pj, j, mimetype)

Print coefficient pj of monomial pj * x^j with the given mimetype.

For pretty printing different number types, or for adding parentheses,
methods can be added to this function. If no mimetype is desired,
adding a method to `Base.show_unquoted` is suggested, as this will
also be useful for the default `show` methods. The following example
shows how `Dual` objects of `DualNumbers` may be printed with
parentheses.

```
using DualNumbers
julia> Poly([Dual(1,2), Dual(3,4)])
Poly(1 + 2ɛ + 3 + 4ɛ*x)
julia> function Base.show_unquoted(io::IO, pj::Dual, indent::Int, prec::Int)
       if Base.operator_precedence(:+) <= prec
               print(io, "(")
               show(io, pj)
               print(io, ")")
           else
               show(io, pj)
           end
end

julia> Poly([Dual(1,2), Dual(3,4)])
Poly((1 + 2ɛ) + (3 + 4ɛ)*x)
```
"""
printcoefficient(io::IO, pj::Any, j, mimetype) = Base.show_unquoted(io, pj, 0, Base.operator_precedence(:*))

# pretty print rational numbers in latex
function printcoefficient(io::IO, a::Rational{T}, j, mimetype::MIME"text/latex") where {T}
    abs(a.den) == one(T) ? print(io, a.num) : print(io, "\\frac{$(a.num)}{$(a.den)}")
end

# print complex numbers with parentheses as needed
function printcoefficient(io::IO, pj::Complex{T}, j, mimetype) where {T}

    hasreal = abs(real(pj)) > 0 || isnan(real(pj)) || isinf(real(pj))
    hasimag = abs(imag(pj)) > 0 || isnan(imag(pj)) || isinf(imag(pj))

    if hasreal && hasimag
        Base.show_unquoted(io, pj, 0, Base.operator_precedence(:*))
    elseif hasreal
        a = real(pj)
        (j==0 || showone(T) || a != one(T)) && printcoefficient(io, a, j, mimetype)
    elseif hasimag
        b = imag(pj)
        (showone(T) || b != one(T)) && printcoefficient(io, b, j,  mimetype)
        (isnan(imag(pj)) || isinf(imag(pj))) && print(io, showop(mimetype, "*"))
        print(io, im)
    else
        return
    end
end



## show exponent

exponent_text(i, ::MIME) = "^$(i)"
exponent_text(i, ::MIME"text/html") = "<sup>$(i)</sup>"
exponent_text(i, ::MIME"text/latex") = "^{$(i)}"

function printexponent(io, var, i, mimetype::MIME)
    if i == 0
        return
    elseif i == 1
        print(io,var)
    else
        print(io, var, exponent_text(i, mimetype))
    end
end


####

## text/plain
Base.show(io::IO, p::Poly{T}) where {T} = show(io, MIME("text/plain"), p)

function Base.show(io::IO, mimetype::MIME"text/plain", p::Poly{T}) where {T}
    print(io,"Poly(")
    printpoly(io, p, mimetype)
    print(io,")")

end

## text/latex
function Base.show(io::IO, mimetype::MIME"text/latex", p::Poly{T}) where {T}
    print(io, "\$")
    printpoly(io, p, mimetype)
    print(io, "\$")
end


## text/html
function Base.show(io::IO, mimetype::MIME"text/html", p::Poly{T}) where {T}
    printpoly(io, p, mimetype)
end
