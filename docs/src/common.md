# Usage

All polynomials have the following functionality. In some cases, there is not a direct function call and therefore the polynomials have to be converted to the standard [`Polynomial`](@ref) type before continuing. 

## Functions

```@index
Pages = ["common.md"]
```

```@docs
AbstractPolynomial
coeffs
order
degree
length
domain
chop
chop!
truncate
truncate!
```

```@docs
+
-
*
/
div
rem
gcd
```

```@docs
fromroots
roots
derivative
integral
integrate
fit
companion
vander
```