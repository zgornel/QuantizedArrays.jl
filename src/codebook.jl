struct CodeBook{U<:Unsigned, T}
    codes::Vector{U}
    vectors::Matrix{T}
    codemap::Dict{U, Int}
end

# Basic parametric constructor
CodeBook(codes::Vector{U}, vectors::Matrix{T}) where {U<:Unsigned, T} = begin
    length(codes) == size(vectors, 2) ||
        throw(DimensionMismatch("Dimension of codes and vectors are inconsistent."))
    codemap = Dict{U, Int}(c => i for (i,c) in enumerate(codes))
    return CodeBook(codes, vectors, codemap)
end

# Show method
Base.show(io::IO, cb::CodeBook{U,T}) where {U,T} = begin
    len_vecs, num_codes = size(cb.vectors)
    print(io, "CodeBook $(num_codes) codes, $(len_vecs)-element $(T) vectors")
end


# Basic indexing: for a key, get the vector
Base.getindex(cb::CodeBook{U,T}, key::U) where {U,T} = cb.vectors[:, cb.codemap[key]]

# Codebook building methods
function build_codebook(aa::AbstractMatrix{T},
                        k::Int,
                        code_type::Type{<:Unsigned};
                        method::Symbol=DEFAULT_METHOD
                       ) where {T}
    if method == :sampling
        vectors = Matrix{T}(undef, size(aa,1), k)
        codes = Vector{code_type}(undef, k)
        cs = floor(Int, size(aa,2)/k)  # column step
        @inbounds @simd for i in 1:k
            cr = cs*(i-1)+1 : cs*i # column range
            codes[i] = code_type(i-1)
            vectors[:,i] = aa[:, rand(cr)]  # randomly sample column
        end
        return CodeBook(codes, vectors)
    elseif method == :pq
        @error "Codebook generation method '$(method)' is not implemented yet."
        # ...
    else
        @error "Unknown codebook generation method '$(method)'"
    end
end


function infer_codes(cb::CodeBook{U,T},
                     aa::AbstractMatrix{T};
                     distance=Distances.Euclidean()
                    ) where {U,T}
    dists = pairwise(distance, cb.vectors, aa, dims=2)  # distances between codebook vectors and data
    best_idx = getindex(findmin(dists, dims=1), 2)      # get positions vector of all minima 
    best_idx_cols = vec(getindex.(best_idx, 1))         # get row values from positions
    return cb.codes[best_idx_cols]                      # get codes corresponding to minima
end
