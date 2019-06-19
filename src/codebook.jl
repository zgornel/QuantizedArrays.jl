struct CodeBook{U<:Unsigned,T}
    codes::Vector{U}
    vectors::Matrix{T}
    codemap::Dict{U,Int}
end


# Constructors
CodeBook(codes::Vector{U}, vectors::Matrix{T}) where {U,T} = begin
    length(codes) != size(vectors, 2) &&
        throw(DimensionMismatch("Dimension of codes and vectors are inconsistent"))
    codemap = Dict{U, Int}(c => i for (i,c) in enumerate(codes))
    return CodeBook(codes, vectors, codemap)
end

CodeBook{U}(vectors::Matrix{T}) where {U,T} = begin
    k = size(vectors,2)
    nb = TYPE_TO_BITS[U]
    nbk = TYPE_TO_BITS[quantized_eltype(k)]
    nbk > nb && throw(ErrorException("$k vectors cannot be coded in $nb bits."))
    codes = Vector{U}(0:k-1)
    codemap = Dict{U, Int}(c => i for (i,c) in enumerate(codes))
    return CodeBook(codes, vectors, codemap)
end

CodeBook{U,T}(vectors::Matrix{T2}) where {U,T<:Integer,T2} =
    CodeBook{U}(round.(T, vectors))

CodeBook{U,T}(vectors::Matrix{T2}) where {U,T<:AbstractFloat,T2} =
    CodeBook{U}(convert.(T, vectors))

CodeBook(vectors::Matrix{T}) where {T} =
    CodeBook{quantized_eltype(size(vectors, 2)), T}(vectors)


# Show method
Base.show(io::IO, cb::CodeBook{U,T}) where {U,T} = begin
    len_vecs, num_codes = size(cb.vectors)
    print(io, "CodeBook{$U,$T} $(num_codes) codes, $(len_vecs)-element vectors")
end


# Basic indexing: for a key, get the vector
Base.getindex(cb::CodeBook{U,T}, key::U) where {U,T} =
    cb.vectors[:, cb.codemap[key]]


# Get compressed codes for an input matrix
function encode(cb::CodeBook{U,T},
                aa::AbstractMatrix{T};
                distance::Distances.PreMetric=DEFAULT_DISTANCE
               ) where {U,T}
    dists = pairwise(distance, cb.vectors, aa, dims=2)  # distances between codebook vectors and data
    best_idx = getindex(findmin(dists, dims=1), 2)      # get positions vector of all minima
    best_idx_cols = vec(getindex.(best_idx, 1))         # get row values from positions
    return cb.codes[best_idx_cols]                      # get codes corresponding to minima
end


# Determine code type
quantized_eltype(k::Int) = begin
    b = clamp(log2(k), 1, MAX_BITS)
    minbits = 1024
    local mintype
    for (nbits, utype) in BITS_TO_TYPE
        if b <= nbits && nbits < minbits
            minbits = nbits
            mintype = utype
        end
    end
    if minbits > MAX_BITS
        @error "Number is too large to fit in $MAX_BITS bits."
    end
    return mintype
end

# Codebooks generation entrypoint
function build_codebooks(X::AbstractMatrix{T},
                         k::Int,
                         m::Int;
                         method::Symbol=DEFAULT_METHOD,
                         distance::Distances.PreMetric=DEFAULT_DISTANCE,
                         kwargs...) where {T}
    nrows, ncols = size(X)
    k = min(k, ncols)                # number of codes for a quantizer
    m = min(m, nrows)                # number of quantizers
    U = quantized_eltype(k)          # type of codes
    if method == :sample
        vectors = sampling_codebooks(X, k, m)  # does not use distances
    elseif method == :pq
        vectors = pq_codebooks(X, k, m, distance=distance; kwargs...)
    elseif method == :opq
        vectors = opq_codebooks(X, k, m, distance=distance; kwargs...)
    else
        @error "Unknown codebook generation method '$(method)'"
    end
    return map(CodeBook{U,T}, vectors)
end
