

macro registerN(name, params...)
    poly = esc(name)
    αs = tuple(esc.(params)...)
    quote
        Base.convert(::Type{P}, q::Q) where {$(αs...),T, P<:$poly{$(αs...),T}, Q <: $poly{$(αs...),T}} = q
        Base.convert(::Type{$poly{$(αs...)}}, q::Q) where {$(αs...),T, Q <: $poly{$(αs...),T}} = q
        Base.promote(p::P, q::Q) where {$(αs...),T, X, P <:$poly{$(αs...),T,X}, Q <: $poly{$(αs...),T,X}} = p,q
        Base.promote_rule(::Type{<:$poly{$(αs...),T,X}}, ::Type{<:$poly{$(αs...),S,X}}) where {$(αs...),T,S,X} =
            $poly{$(αs...),promote_type(T, S),X}
        Base.promote_rule(::Type{<:$poly{$(αs...),T,X}}, ::Type{S}) where {$(αs...),T,X,S<:Number} =
            $poly{$(αs...),promote_type(T,S),X}

        function $poly{$(αs...),T}(x::AbstractVector{S}, var::SymbolLike = Var(:x)) where {$(αs...),T,S}
            $poly{$(αs...),T,Symbol(var)}(T.(x))
        end
        $poly{$(αs...)}(coeffs::AbstractVector{T}, var::SymbolLike=Var(:x)) where {$(αs...),T} =
            $poly{$(αs...),T,Symbol(var)}(coeffs)
        $poly{$(αs...),T,X}(c::AbstractPolynomial{S,Y}) where {$(αs...),T,X,S,Y} = convert($poly{$(αs...),T,X}, c)
        $poly{$(αs...),T}(c::AbstractPolynomial{S,Y}) where {$(αs...),T,S,Y} = convert($poly{$(αs...),T}, c)
        $poly{$(αs...),}(c::AbstractPolynomial{S,Y}) where {$(αs...),S,Y} = convert($poly{$(αs...),}, c)

        $poly{$(αs...),T,X}(n::Number) where {$(αs...),T,X} = T(n)*one($poly{$(αs...),T,X})
        $poly{$(αs...),T}(n::Number, var::SymbolLike = Var(:x)) where {$(αs...),T} = T(n)*one($poly{$(αs...),T,Symbol(var)})
        $poly{$(αs...)}(n::S, var::SymbolLike = Var(:x)) where {$(αs...), S<:Number} =
            n*one($poly{$(αs...),S,Symbol(var)})
        $poly{$(αs...),T}(var::SymbolLike=Var(:x)) where {$(αs...), T} =
            variable($poly{$(αs...),T,Symbol(var)})
        $poly{$(αs...)}(var::SymbolLike=Var(:x)) where {$(αs...)} = variable($poly{$(αs...)},Symbol(var))
        (p::$poly)(x) = evalpoly(x, p)
    end
end
