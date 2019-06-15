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
function infer_codes(cb::CodeBook{U,D,T}, aa::AbstractMatrix{T}) where {U,D,T}
    dists = pairwise(cb.distance, cb.vectors, aa, dims=2)  # distances between codebook vectors and data
    best_idx = getindex(findmin(dists, dims=1), 2)         # get positions vector of all minima 
    best_idx_cols = vec(getindex.(best_idx, 1))            # get row values from positions
    return cb.codes[best_idx_cols]                         # get codes corresponding to minima
end


# Codebook building methods
function build_codebook(aa::AbstractMatrix{T},
                        k::Int,
                        code_type::Type{<:Unsigned};
                        method::Symbol=DEFAULT_METHOD,
                        distance::Distances.PreMetric=DEFAULT_DISTANCE
                       ) where {T}
    if method == :sample
        return _build_sampling_codebook(aa, k, code_type)
    elseif method == :pq
        return _build_kmeans_codebook(aa, k, code_type, distance=distance)
    else
        @error "Unknown codebook generation method '$(method)'"
    end
end


function _build_sampling_codebook(aa::AbstractMatrix{T},
                                  k::Int,
                                  code_type::Type{<:Unsigned};
                                 ) where {T}
    vectors = Matrix{T}(undef, size(aa,1), k)
    codes = Vector{code_type}(undef, k)
    cs = floor(Int, size(aa,2)/k)  # column step
    @inbounds @simd for i in 1:k
        cr = cs*(i-1)+1 : cs*i # column range
        codes[i] = code_type(i-1)
        vectors[:,i] = aa[:, rand(cr)]  # randomly sample column
    end
    return CodeBook(codes, vectors)
end

function _build_kmeans_codebook(aa::AbstractMatrix{T},
                                k::Int,
                                code_type::Type{<:Unsigned};
                                distance::Distances.PreMetric=DEFAULT_DISTANCE
                               ) where {T}
    codes = Vector{code_type}(0:k-1)
    model = kmeans(aa, k, maxiter=30, display=:none)
    return CodeBook(codes, T.(model.centers), distance=distance)
end
