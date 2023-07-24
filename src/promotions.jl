
# do both ways to avoid an issue with the next set of promotion rules
Base.promote_rule(::Type{P}, ::Type{Q}) where {B,T,S,X,
                                               P<:MutableDensePolynomial{B,T,X},
                                               Q<:MutableSparsePolynomial{B,S,X}} = MutableDensePolynomial{B,promote_type(T,S),X}
Base.promote_rule(::Type{Q}, ::Type{P}) where {B,T,S,X,
                                               P<:MutableDensePolynomial{B,T,X},
                                               Q<:MutableSparsePolynomial{B,S,X}} = MutableDensePolynomial{B,promote_type(T,S),X}
Base.promote_rule(::Type{P}, ::Type{Q}) where {B,T,S,X,N,
                                               P<:MutableDensePolynomial{B,T,X},
                                               Q<:ImmutableDensePolynomial{B,S,X,N}} = MutableDensePolynomial{B,promote_type(T,S),X}
Base.promote_rule(::Type{Q}, ::Type{P}) where {B,T,S,X,N,
                                               P<:MutableDensePolynomial{B,T,X},
                                               Q<:ImmutableDensePolynomial{B,S,X,N}} = MutableDensePolynomial{B,promote_type(T,S),X}
Base.promote_rule(::Type{P}, ::Type{Q}) where {B,T,S,X,N,
                                               P<:MutableSparsePolynomial{B,T,X},
                                               Q<:ImmutableDensePolynomial{B,S,X,N}} = MutableSparsePolynomial{B,promote_type(T,S),X}
Base.promote_rule(::Type{Q}, ::Type{P}) where {B,T,S,X,N,
                                               P<:MutableSparsePolynomial{B,T,X},
                                               Q<:ImmutableDensePolynomial{B,S,X,N}} = MutableSparsePolynomial{B,promote_type(T,S),X}


# XXX need both to work around more general promotion to Polynomial type
# Base.promote_rule(::Type{P},::Type{Q}) where {B<:StandardBasis,T,X, P<:AbstractUnivariatePolynomial{B,T,X},
#                                               S,   Q<:AbstractPolynomial{S,X}} =
#                                                   MutableDensePolynomial{StandardBasis, promote_type(T, S), X}
# Base.promote_rule(::Type{Q},::Type{P}) where {B<:StandardBasis,T,X, P<:AbstractUnivariatePolynomial{B,T,X},
#                                               S,   Q<:AbstractPolynomial{S,X}} =
#                                                   MutableDensePolynomial{StandardBasis, promote_type(T, S), X}


# Methods to ensure that matrices of polynomials behave as desired
Base.promote_rule(::Type{P},::Type{Q}) where {T,X, P<:AbstractPolynomial{T,X},
                                              S,   Q<:AbstractPolynomial{S,X}} =
                                                   Polynomial{promote_type(T, S),X}

Base.promote_rule(::Type{P},::Type{Q}) where {T,X, P<:AbstractPolynomial{T,X},
                                              S,Y, Q<:AbstractPolynomial{S,Y}} =
                                                  assert_same_variable(X,Y)
