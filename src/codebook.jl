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


# Codebook building algorithms
# ----------------------------
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
        return build_sample_codebooks(X, k, m, U)  # does not use distances
    elseif method == :pq
        return build_pq_codebooks(X, k, m, U, distance=distance; kwargs...)
    elseif method == :opq
        return build_opq_codebooks(X, k, m, U, distance=distance; kwargs...)
    else
        @error "Unknown codebook generation method '$(method)'"
    end
end


# Utility function that returns a row ranged based on an iteration index
@inline rowrange(nrows::Int, m::Int, i::Int) = begin
    rs = floor(Int, nrows/m)  # row step
    rr = rs*(i-1)+1 : rs*i    # row range
    return rr
end


function build_sample_codebooks(X::AbstractMatrix{T},
                                k::Int,
                                m::Int,
                                ::Type{U}) where {U,T}
    nrows, ncols = size(X)
    cbooks = Vector{CodeBook{U,T}}(undef, m)
    @inbounds for i in 1:m
        rr = rowrange(nrows, m, i)
        vectors = Matrix{T}(X[rr, sample(1:ncols, k, replace=false)])
        cbooks[i] = CodeBook{U,T}(vectors)
    end
    return cbooks
end

function build_pq_codebooks(X::AbstractMatrix{T},
                            k::Int,
                            m::Int,
                            ::Type{U};
                            distance::Distances.PreMetric=DEFAULT_DISTANCE,
                            maxiter::Int=DEFAULT_PQ_MAXITER
                           ) where {U,T}
    nrows = size(X, 1)
    cbooks = Vector{CodeBook{U,T}}(undef, m)
    @inbounds for i in 1:m
        rr = rowrange(nrows, m, i)
        model = kmeans(X[rr, :], k,
                       maxiter=maxiter,
                       distance=distance,
                       init=:kmpp,
                       display=:none)
        cbooks[i] = CodeBook{U,T}(model.centers)
    end
    return cbooks
end

function build_opq_codebooks(X::AbstractMatrix{T},
                             k::Int,
                             m::Int,
                             ::Type{U};
                             distance::Distances.PreMetric=DEFAULT_DISTANCE,
                             maxiter::Int=DEFAULT_OPQ_MAXITER
                            ) where {U,T}
    # Initialize R
    nrows, ncols = size(X)
    R = diagm(0 => ones(T, nrows))
    X̂ = R' * X

    # Initialize codebooks, codes
    codes = fill(zeros(Int, ncols), m)
    vectors = [X̂[rowrange(nrows, m, i), sample(1:ncols, k, replace=false)]
               for i in 1:m]
    cweights   = zeros(Int, k)
    costs      = zeros(ncols)
    counts     = zeros(Int, k)
    unused     = Int[]
    to_update = zeros(Bool, k)
    @inbounds for i in 1:m
        rr = rowrange(nrows, m, i)
        dists = Distances.pairwise(distance, vectors[i], X̂[rr, :], dims=2)
        Clustering.update_assignments!(dists, false, codes[i], costs, counts, to_update, unused)
        X̂[rr,:] .= vectors[i][:, codes[i]]
    end

    # Run optimization
    to_update = ones(Bool, k)
    for _ = 1:maxiter
        # Update R using orthogonal Procustes closed form solution
        # and update data rotated data matrix X̂
        Uₓ,_,Vₓ = svd(X * X̂', full=false)
        R = Vₓ * Uₓ'
        X̂ = R * X
        @inbounds for i in 1:m
            rr = rowrange(nrows, m, i)
            # Update subspace cluster centers
            Clustering.update_centers!(X̂[rr, :], nothing, codes[i], to_update, vectors[i], cweights)
            # Update subspace data assignments
            dists = Distances.pairwise(distance, vectors[i], X̂[rr, :], dims=2)
            Clustering.update_assignments!(dists, false, codes[i], costs, counts, to_update, unused)
            X̂[rr, :] .= vectors[i][:, codes[i]]
        end
    end
    return map(CodeBook{U,T}, vectors)
end
