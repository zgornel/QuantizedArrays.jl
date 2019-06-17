struct CodeBook{U<:Unsigned,T}
    codes::Vector{U}
    vectors::Matrix{T}
    codemap::Dict{U,Int}
end

# Basic parametric constructor
CodeBook(codes::Vector{U}, vectors::Matrix{T}) where {U,T} = begin
    length(codes) == size(vectors, 2) ||
        throw(DimensionMismatch("Dimension of codes and vectors are inconsistent"))
    codemap = Dict{U, Int}(c => i for (i,c) in enumerate(codes))
    return CodeBook(codes, vectors, codemap)
end


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


# Codebook building methods
function build_codebooks(aa::AbstractMatrix{T},
                         k::Int,
                         m::Int;
                         method::Symbol=DEFAULT_METHOD,
                         distance::Distances.PreMetric=DEFAULT_DISTANCE
                        ) where {T}
    nrows, ncols = size(aa)
    k = min(k, ncols)                # number of codes for a quantizer
    m = min(m, nrows)                # number of quantizers
    U = quantized_eltype(k)          # type of codes

    if method == :sample
        return build_sample_codebooks(aa, k, m, U)  # does not use distances
    elseif method == :pq
        return build_pq_codebooks(aa, k, m, U, distance=distance)
    elseif method == :opq
        return build_opq_codebooks(aa, k, m, U, distance=distance)
    else
        @error "Unknown codebook generation method '$(method)'"
    end
end


function build_sample_codebooks(aa::AbstractMatrix{T},
                                k::Int,
                                m::Int,
                                ctype::Type{U}) where {U,T}
    nrows, ncols = size(aa)
    cs = floor(Int, ncols/k)  # column step
    rs = floor(Int, nrows/m)  # row step

    cbooks = Vector{CodeBook{U,T}}(undef, m)
    @inbounds @simd for i in 1:m
        rr = rs*(i-1)+1 : rs*i  # row range
        codes = Vector{U}(0:k-1)
        vectors = Matrix{T}(aa[rr, sample(1:ncols, k, replace=false)])
        cbooks[i] = CodeBook(codes, vectors)
    end
    return cbooks
end

function build_pq_codebooks(aa::AbstractMatrix{T},
                            k::Int,
                            m::Int,
                            ctype::Type{U};
                            distance::Distances.PreMetric=DEFAULT_DISTANCE
                           ) where {U,T}
    nrows, ncols = size(aa)
    cs = floor(Int, ncols/k)  # column step
    rs = floor(Int, nrows/m)  # row step

    cbooks = Vector{CodeBook{U,T}}(undef, m)
    @inbounds @simd for i in 1:m
        rr = rs*(i-1)+1 : rs*i  # row range
        codes = Vector{U}(0:k-1)
        model = kmeans(aa[rr, :], k, maxiter=30, init=:kmpp, display=:none)
        cbooks[i] = CodeBook(codes, T.(model.centers))
    end
    return cbooks
end

function build_opq_codebooks(aa::AbstractMatrix{T},
                             k::Int,
                             m::Int,
                             ctype::Type{U};
                             distance::Distances.PreMetric=DEFAULT_DISTANCE
                            ) where {U,T}
    nrows, ncols = size(aa)
    cs = floor(Int, ncols/k)  # column step
    rs = floor(Int, nrows/m)  # row step

    cbooks = Vector{CodeBook{U,T}}(undef, m)
    # Initialize R, codebooks, codes
    # repeat
    #   Step(1): project the data X̂ = RX
    #   for i = 1:m
    #       for j = 1:k update ĉᵐ(j) by the sample mean of x̂ᵐ i.e. subspace cluster centers by projected samples form the subspace
    #       for any x̂ᵐ update iᵐ(x̂ᵐ) i.e. codes by the subindex of the codeword ĉᵐ that is nearest to x̂ᵐ
    #   end
    #   Solve R by
    #     U,S,V = svd(X*Xq', full=false)
    #     R = V*U'
    #until max number of iterations
    return cbooks
end
