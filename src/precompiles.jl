# precompiles

let p = fromroots(Polynomial, [1,1,2])
    Multroot.multroot(p)
    gcd(p, derivative(p); method=:numerical)
    #uvw(p, derivative(p); method=:numerical)
end
