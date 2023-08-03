Base.promote_rule(::Type{P}, ::Type{Q}) where {B,T,S,X,
                                               P<:AbstractUnivariatePolynomial{B,T,X},
                                               Q<:AbstractUnivariatePolynomial{B,S,X}} = MutableDenseLaurentPolynomial{B,promote_type(T,S),X}

Base.promote_rule(::Type{P}, ::Type{Q}) where {B,T,S,X,
                                               P<:AbstractPolynomial{T,X},
                                               Q<:AbstractUnivariatePolynomial{B,S,X}} = MutableDenseLaurentPolynomial{B,promote_type(T,S),X}

Base.promote_rule(::Type{P}, ::Type{Q}) where {B,T,S,X,
                                               P<:AbstractUnivariatePolynomial{B,T,X},
                                               Q<:AbstractPolynomial{S,X}} = MutableDenseLaurentPolynomial{B,promote_type(T,S),X}

## XXX these are needed for rational-functions
# Base.promote_rule(::Type{P}, ::Type{Q}) where {B,T,S,X,
#                                                P<:Polynomial{T,X},
#                                                Q<:AbstractUnivariatePolynomial{B,S,X}} = Polynomial{promote_type(T,S),X}

# Base.promote_rule(::Type{P}, ::Type{Q}) where {B,T,S,X,
#                                                P<:AbstractUnivariatePolynomial{B,T,X},
#                                                Q<:Polynomial{S,X}} = Polynomial{promote_type(T,S),X}

# # L,P -> L (also S,P -> L should be defined...)
# Base.promote_rule(::Type{P}, ::Type{Q}) where {T,S,X,
#                                                P<:Polynomial{T,X},
#                                                Q<:LaurentPolynomial{S,X}} = LaurentPolynomial{promote_type(T,S),X}

# Base.promote_rule(::Type{P}, ::Type{Q}) where {T,S,X,
#                                                P<:LaurentPolynomial{T,X},
#                                                Q<:Polynomial{S,X}} = LaurentPolynomial{promote_type(T,S),X}



# Methods to ensure that matrices of polynomials behave as desired
Base.promote_rule(::Type{P},::Type{Q}) where {T,X, P<:AbstractPolynomial{T,X},
                                              S,   Q<:AbstractPolynomial{S,X}} =
                                                   Polynomial{promote_type(T, S),X}

Base.promote_rule(::Type{P},::Type{Q}) where {T,X, P<:AbstractPolynomial{T,X},
                                              S,Y, Q<:AbstractPolynomial{S,Y}} =
                                                  assert_same_variable(X,Y)
