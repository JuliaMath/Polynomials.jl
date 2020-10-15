export printpoly

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
hasneg(::Type{<:AbstractPolynomial{S}}) where {S} = false
showone(::Type{<:AbstractPolynomial{S}}) where {S} = false



#=
Common Printing
=#

Base.show(io::IO, p::AbstractPolynomial) = show(io, MIME("text/plain"), p)

function Base.show(io::IO, mimetype::MIME"text/plain", p::P) where {P<:AbstractPolynomial}
    print(io,"$(P.name)(")
    printpoly(IOContext(io, :compact=>get(io, :compact, false)), p, mimetype)
    print(io,")")
end

function Base.show(io::IO, mimetype::MIME"text/latex", p::AbstractPolynomial)
    print(io, "\$")
    printpoly(io, p, mimetype)
    print(io, "\$")
end

function Base.show(io::IO, mimetype::MIME"text/html", p::AbstractPolynomial)
    printpoly(io, p, mimetype)
end

#####

"Show different operations depending on mimetype. `l-` is leading minus sign."
function showop(::MIME"text/plain", op)
    d = Dict("*" => "*", "+" => " + ", "-" => " - ", "l-" => "-")
    get(d, op, "")
end

function showop(::MIME"text/latex", op)
    d = Dict("*" => "\\cdot ", "+" => " + ", "-" => " - ", "l-" => "-")
    get(d, op, "")
end

function showop(::MIME"text/html", op)
    d = Dict("*" => "&#8729;", "+" => " &#43; ", "-" => " &#45; ", "l-" => "&#45;")
    get(d, op, "")
end



###

"""
    printpoly(io::IO, p::AbstractPolynomial, mimetype = MIME"text/plain"(); descending_powers=false, offset::Int=0)

Print a human-readable representation of the polynomial `p` to `io`. The MIME
types "text/plain" (default), "text/latex", and "text/html" are supported. By
default, the terms are in order of ascending powers, matching the order in
`coeffs(p)`; specifying `descending_powers=true` reverses the order.
`offset` allows for an integer number to be added to the exponent, just for printing.
`var` allows for overriding the variable used for printing. Setting multiplication symbol to `""`
will avoid an operator being printed. Setting `compact=true` will use a compact style for floating point numbers.

# Examples
```jldoctest show
julia> using Polynomials

julia> printpoly(stdout, Polynomial([1,2,3], :y))
1 + 2*y + 3*y^2

julia> printpoly(stdout, Polynomial([1,2,3], :y), descending_powers=true)
3*y^2 + 2*y + 1

julia> printpoly(stdout, Polynomial([2, 3, 1], :z), descending_powers=true, offset=-2)
1 + 3*z^-1 + 2*z^-2

julia> printpoly(stdout, Polynomial([-1, 0, 1], :z), offset=-1, descending_powers=true)
z - z^-1

julia> printpoly(stdout, Polynomial([-1, 0, 1], :z), offset=-1, descending_powers=true, var=:x)
x - x^-1
```
"""
function printpoly(io::IO, p::P, mimetype=MIME"text/plain"();
                   descending_powers=false, offset::Int=0, var=p.var,
                   compact=false, mulsymbol="*") where {T,P<:AbstractPolynomial{T}}
    first = true
    printed_anything = false
    for i in (descending_powers ? reverse(eachindex(p)) : eachindex(p))
        ioc = IOContext(io, :compact=>get(io, :compact, compact),
                        :multiplication_symbol => get(io, :multiplication_symbol, mulsymbol))
        printed = showterm(ioc, P, p[i], var, i+offset, first, mimetype)
        first &= !printed
        printed_anything |= printed
    end
    printed_anything || print(io, zero(T))
    return nothing
end

"""
    showterm(io::IO, ::Type{<:AbstractPolynomial} pj, var, j, first, mimetype)

Shows the j'th term of the given polynomial. Returns `true` after successfully printing.

For example. for a `Polynomial` this would show the term `pj * var^j`.
"""
function showterm(io::IO, ::Type{AbstractPolynomial}, pj::T, var, j, first::Bool, mimetype) where {T} end


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
## pass `:multiplication_symbol => "" to IOContext have no sign
function printproductsign(io::IO, pj::T, j, mimetype) where {T}
    j == 0 && return
    multiplication_symbol = showop(mimetype, get(io, :multiplication_symbol,"*"))
    (showone(T) || pj != one(T)) &&  print(io, multiplication_symbol)
end

function printproductsign(io::IO, pj::T, j, mimetype) where {T<:Complex}
    j == 0 && return
    (a,b) = reim(pj)
    !iszero(a) && !iszero(b) && return # parentheses inserted, no * needed
    !iszero(a) && return printproductsign(io, a, j, mimetype)
    print(io, showop(mimetype, "*"))
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
julia> Polynomial([Dual(1,2), Dual(3,4)])
Polynomial(1 + 2ɛ + 3 + 4ɛ*x)

julia> function Base.show_unquoted(io::IO, pj::Dual, indent::Int, prec::Int)
       if Base.operator_precedence(:+) <= prec
            print(io, "(")
            show(io, pj)
            print(io, ")")
        else
            show(io, pj)
        end
    end

julia> Polynomial([Dual(1,2), Dual(3,4)])
Polynomial((1 + 2ɛ) + (3 + 4ɛ)*x)
```
"""
printcoefficient(io::IO, pj::Any, j, mimetype) = Base.show_unquoted(io, pj, 0, Base.operator_precedence(:*))

# pretty print rational numbers in latex
function printcoefficient(io::IO, a::Rational{T}, j, mimetype::MIME"text/latex") where {T}
    abs(a.den) == one(T) ? print(io, a.num) : print(io, "\\frac{$(a.num)}{$(a.den)}")
end

# print complex numbers with parentheses as needed
function printcoefficient(io::IO, pj::S, j, mimetype) where {T,S <: Complex{T}}

    (a,b) = reim(pj)
    hasreal = !iszero(a) || isnan(a) || isinf(a)
    hasimag = !iszero(b) || isnan(b) || isinf(b)

    if hasreal && hasimag
        iszero(j) || print(io, "(")
        print(io, a)

        # print b
        if isone(b) || isone(-b)
            if hasneg(S) && b < 0
                print(io, showop(mimetype, "-"))
            else
                print(io, showop(mimetype, "+"))
            end
        else
            if hasneg(S) && b < 0
                print(io, showop(mimetype, "-"))
                (showone(S) || !isone(-b)) && print(io, -b)
            else
                print(io, showop(mimetype,"+"))
                print(io, b)
            end
            (isnan(b) || isinf(b)) && print(io, showop(mimetype, "*"))
        end

        print(io, imagsymbol(mimetype))
        iszero(j) || print(io, ")")

    elseif hasreal

        (iszero(j) || showone(T) || !isone(a)) && printcoefficient(io, a, j, mimetype)

    elseif hasimag
        (showone(T) || !isone(b)) && printcoefficient(io, b, j, mimetype)
        (isnan(b) || isinf(b)) && print(io, showop(mimetype, "*"))
        print(io, imagsymbol(mimetype))

    end

    return nothing

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

function unicode_exponent(io, j)
    a = ("⁻","","","⁰","¹","²","³","⁴","⁵","⁶","⁷","⁸","⁹")
    for i in string(j)
        print(io, a[Int(i)-44])
    end
end

function unicode_subscript(io, j)
    a = ("₋","","","₀","₁","₂","₃","₄","₅","₆","₇","₈","₉")
    for i in string(j)
        print(io, a[Int(i)-44])
    end
end

imagsymbol(::Any) = "im"
imagsymbol(::MIME"text/latex") = "i"
