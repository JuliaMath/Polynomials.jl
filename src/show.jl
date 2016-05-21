
function show(io::IO, p::Poly)
    print(io,"Poly(")
    print(io,p)
    print(io,")")
end

function writemime{T<:Number}(io::IO, ::MIME"text/latex",
                              poly::Polynomials.Poly{T})
    print(io, "\$")
    print(io, poly, true)
    print(io, "\$")
end

function writemime{T<:Number}(io::IO, ::MIME"text/latex",
                              pade::Polynomials.Pade{T,T})
    print(io, "\$\$\\frac{")
    print(io, pade.p, true)
    print(io, "}{")
    print(io, pade.q, true)
    print(io, "}\$\$")
end

function printexponent(io,var,i)
    if i == 0
    elseif i == 1
        print(io,var)
    else
        print(io,var,"^",i)
    end
end

function show(io::IO, r::Rational, latex_print::Bool)
    if latex_print
        if den(r) == 1
            print(io, string(num(r)))
        else
            print(io, "\\frac{", string(num(r)), "}{", string(den(r)), "}")
        end
    else
        print(io, string(r))
    end
end
show{N<:Number}(io::IO, n::N, latex_print::Bool) = show(io, n)

function printterm{T}(io::IO,p::Poly{T},j,first,latex_print::Bool)
    pj = p[j]
    if pj == zero(T)
        return false
    end
    neg = pj < zero(T)
    if first
        neg && print(io, "-")    #Prepend - if first and negative
    else
        neg ? print(io, " - ") : print(io," + ")
    end
    pj = abs(pj)
    if pj != one(T) || j == 0
        show(io,pj,latex_print)
    end
    printexponent(io,p.var,j)
    true
end

function printterm{T<:Complex}(io::IO,p::Poly{T},j,first,latex_print::Bool)
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
            print(io,'(', abs(imag(pj)),"im)")
        end
    end
    printexponent(io,p.var,j)
    true
end

function print{T}(io::IO, p::Poly{T}, latex_print::Bool = false)
    first = true
    printed_anything = false
    n = length(p)-1
    for i = 0:n
        printed = printterm(io,p,i,first,latex_print)
        first &= !printed
        printed_anything |= printed
    end
    printed_anything || print(io,zero(T))
end
