struct CodeBook{U<:Unsigned,D<:Distances.PreMetric,T}
    codes::Vector{U}
    vectors::Matrix{T}
    codemap::Dict{U,Int}
    distance::D
end

# Basic parametric constructor
CodeBook(codes::Vector{U},
         vectors::Matrix{T};
         distance::Distances.PreMetric=DEFAULT_DISTANCE
        ) where {U<:Unsigned, T} = begin
    length(codes) == size(vectors, 2) ||
        throw(DimensionMismatch("Dimension of codes and vectors are inconsistent."))
    codemap = Dict{U, Int}(c => i for (i,c) in enumerate(codes))
    return CodeBook(codes, vectors, codemap, distance)
end


# Show method
Base.show(io::IO, cb::CodeBook{U,D,T}) where {U,D,T} = begin
    len_vecs, num_codes = size(cb.vectors)
    print(io, "CodeBook $(num_codes) codes, $(len_vecs)-element $T vectors, $D distance")
end


# Basic indexing: for a key, get the vector
Base.getindex(cb::CodeBook{U,D,T}, key::U) where {U,D,T} =
    cb.vectors[:, cb.codemap[key]]


# Get compressed codes for an input matrix
function encode(cb::CodeBook{U,D,T}, aa::AbstractMatrix{T}) where {U,D,T}
    dists = pairwise(cb.distance, cb.vectors, aa, dims=2)  # distances between codebook vectors and data
    best_idx = getindex(findmin(dists, dims=1), 2)         # get positions vector of all minima 
    best_idx_cols = vec(getindex.(best_idx, 1))            # get row values from positions
    return cb.codes[best_idx_cols]                         # get codes corresponding to minima
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
    D = typeof(distance)             # type of distance
    U = quantized_eltype(k)          # type of codes
    cbooks = Vector{CodeBook{U,D,T}}(undef, m)

    if method == :sample
        update_sample_codebooks!(cbooks, aa, k)
    elseif method == :pq
        update_pq_codebooks!(cbooks, aa, k)
    elseif method == :opq
        update_opq_codebooks!(cbooks, aa, k)
    else
        @error "Unknown codebook generation method '$(method)'"
    end
    return cbooks
end


function update_sample_codebooks!(codebooks::Vector{CodeBook{U,D,T}},
                                  aa::AbstractMatrix{T},
                                  k::Int) where {U,D,T}
    nrows, ncols = size(aa)
    m = length(codebooks)
    cs = floor(Int, ncols/k)  # column step
    rs = floor(Int, nrows/m)  # row step

    @inbounds @simd for i in 1:m
        rr = rs*(i-1)+1 : rs*i  # row range
        codes = Vector{U}(0:k-1)
        vectors = Matrix{T}(aa[rr, sample(1:ncols, k, replace=false)])
        codebooks[i] = CodeBook(codes, vectors, distance=D())
    end
end

function update_pq_codebooks!(codebooks::Vector{CodeBook{U,D,T}},
                              aa::AbstractMatrix{T},
                              k::Int) where {U,D,T}
    nrows, ncols = size(aa)
    m = length(codebooks)
    cs = floor(Int, ncols/k)  # column step
    rs = floor(Int, nrows/m)  # row step

    @inbounds @simd for i in 1:m
        rr = rs*(i-1)+1 : rs*i  # row range
        codes = Vector{U}(0:k-1)
        model = kmeans(aa[rr, :], k, maxiter=30, init=:kmpp, display=:none)
        codebooks[i] = CodeBook(codes, T.(model.centers), distance=D())
    end
end

function update_opq_codebooks!(codebooks::Vector{CodeBook{U,D,T}},
                               aa::AbstractMatrix{T},
                               k::Int) where {U,D,T}
    #TODO: Implement
    @error "Not implemented"
end
