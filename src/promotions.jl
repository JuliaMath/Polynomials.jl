Base.promote_rule(::Type{P}, ::Type{Q}) where {B,T,S,X,
                                               P<:AbstractUnivariatePolynomial{B,T,X},
                                               Q<:AbstractUnivariatePolynomial{B,S,X}} = MutableDenseLaurentPolynomial{B,promote_type(T,S),X}

Base.promote_rule(::Type{P}, ::Type{Q}) where {B,T,S,X,
                                               P<:AbstractDenseUnivariatePolynomial{B,T,X},
                                               Q<:AbstractDenseUnivariatePolynomial{B,S,X}} = MutableDensePolynomial{B,promote_type(T,S),X}

Base.promote_rule(::Type{P}, ::Type{Q}) where {B,T,S,X,
                                               P<:AbstractDenseUnivariatePolynomial{B,T,X},
                                               Q<:AbstractLaurentUnivariatePolynomial{B,S,X}} = MutableDenseLaurentPolynomial{B,promote_type(T,S),X}

Base.promote_rule(::Type{P}, ::Type{Q}) where {B,T,S,X,
                                               P<:AbstractLaurentUnivariatePolynomial{B,T,X},
                                               Q<:AbstractLaurentUnivariatePolynomial{B,S,X}} = MutableDenseLaurentPolynomial{B,promote_type(T,S),X}


# Methods to ensure that matrices of polynomials behave as desired
Base.promote_rule(::Type{P},::Type{Q}) where {T,X, P<:AbstractPolynomial{T,X},
                                              S,   Q<:AbstractPolynomial{S,X}} =
                                                   LaurentPolynomial{promote_type(T, S),X}


Base.promote_rule(::Type{P},::Type{Q}) where {T,X, P<:AbstractPolynomial{T,X},
                                              S,Y, Q<:AbstractPolynomial{S,Y}} =
                                                  assert_same_variable(X,Y)
