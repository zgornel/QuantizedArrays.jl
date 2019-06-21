"""
    CodeBook{U<:Unsigned,T}

The codebook structure. It holds the codes corresponding
to the vector prototypes and the mapping bethween the codes
and prototypes.

# Fields
  * `codes::Vector{U}` the codes
  * `vectors::Matrix{T}` the prototypes
  * `codemap::Dict{U,Int}` mapping from code to column in `vectors`
"""
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


"""
    encode(codebook, aa [;distance=DEFAULT_DISTANCE])

Encodes the input array aa using `distance` to calculate the
closest vector prototype from the `codebook`.
"""
function encode(cb::CodeBook{U,T},
                aa::AbstractMatrix{T};
                distance::Distances.PreMetric=DEFAULT_DISTANCE
               ) where {U,T}
    dists = pairwise(distance, cb.vectors, aa, dims=2)  # distances between codebook vectors and data
    best_idx = getindex(findmin(dists, dims=1), 2)      # get positions vector of all minima
    best_idx_cols = vec(getindex.(best_idx, 1))         # get row values from positions
    return cb.codes[best_idx_cols]                      # get codes corresponding to minima
end


"""
    build_codebooks(X, k, m, U [;method=DEFAULT_METHOD, distance=DEFAULT_DISTANCE, kwargs])

Generates `m` codebooks of `k` prototypes each for the input matrix `X`
using the algorithm specified my `method` and distance `distance`.
Specific codebook aglorithm keyword arguments can be specified as well.

# Arguments
  * `X::AbstractMatrix{T}` input matrix of type `T`
  * `k::Int` number of prototypes/codebook
  * `m::Int` number of codebooks
  * `U::Type{<:Unsigned}` type for codebook codes

# Keyword arguments
  * `method::Symbol` the algorithm to be employed for codebook
generation; possible values are `:sample` (default), `:pq` for
classical k-means clustering codebooks and `:opq` for 'cartesian'
k-means clustering codebooks
  * `distance::PreMetric` the distance to be used in the
codebook generation methods and data encoding
  * `kwargs...` other codebook generation algorithm specific
keyword arguments such as `maxiter::Int`.
"""
function build_codebooks(X::AbstractMatrix{T},
                         k::Int,
                         m::Int,
                         ::Type{U};
                         method::Symbol=DEFAULT_METHOD,
                         distance::Distances.PreMetric=DEFAULT_DISTANCE,
                         kwargs...) where {U,T}
    nrows, ncols = size(X)
    k = min(k, ncols-1)              # number of codes for a quantizer
    m = min(m, nrows)                # number of quantizers
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
