struct QuantizedArray{U<:Unsigned,D<:Distances.PreMetric,T,N} <: AbstractArray{T,N}
    data::Array{U,N}
    quantizer::ArrayQuantizer{U,D,T,N}
end


# Aliases
const QuantizedVector{U,D,T} = QuantizedArray{U,D,T,1}
const QuantizedMatrix{U,D,T} = QuantizedArray{U,D,T,2}


# Basic constructor
function QuantizedArray(aa::AbstractArray{T,N};
                        k::Int=DEFAULT_K,
                        m::Int=DEFAULT_M,
                        method::Symbol=DEFAULT_METHOD,
                        distance::Distances.PreMetric=DEFAULT_DISTANCE
                       ) where {T,N}
    @assert N <=2 "Array quantization is supported only for Vectors and Matrices"
    @assert k >= 1 "`k` has to be larger or equal to 1"
    @assert m >= 1 "`m` has to be larger or equal to 1"
    if N == 2
        nrows = size(aa, 1)
        @assert rem(nrows, m) == 0 "`m` has to divide exactly $nrows rows"
    end
    quantizer = build_quantizer(aa, k=k, m=m, method=method, distance=distance)
    qa = quantize(quantizer, aa)
    return QuantizedArray(qa, quantizer)
end


# eltype, size, length
Base.eltype(qa::QuantizedArray{U,D,T,N}) where {U,D,T,N} = T
Base.size(qa::QuantizedArray) = qa.quantizer.dims
Base.IndexStyle(::Type{<:QuantizedArray}) = IndexLinear()
quantizer(qa::QuantizedArray) = qa.quantizer


# Indexing interface: getindex, setindex!
@inline function Base.getindex(qa::QuantizedVector, i::Int)
    @boundscheck checkbounds(qa.data, i)
    qkey = getindex(qa.data, i)   # get quantization code
    cb = codebooks(qa.quantizer)  # get codebooks
    return cb[qkey][1]            # get quantized value
end

@inline function Base.getindex(qa::QuantizedMatrix, i::Int)
    cbooks = codebooks(quantizer(qa))
    m = length(cbooks)
    col, row, cidx = indices(size(qa, 1), m, i)  # quantized index, partial index
    qkey = getindex(qa.data, cidx, col)          # get quantization code
    return cbooks[cidx][qkey][row]               # get quantized value
end

@inline function Base.getindex(qa::QuantizedMatrix, i::Int, j::Int)
    cbooks = codebooks(quantizer(qa))
    m = length(cbooks)
    _, row, cidx = indices(size(qa, 1), m, i)  # quantized index, partial index
    qkey = getindex(qa.data, cidx, j)          # get quantization code
    return cbooks[cidx][qkey][row]             # get quantized value
end

function indices(n::Int, m::Int, I::Vararg{Int,N}) where {N}
    step = Int(n/m)
    col, row = divrem(I[1]-1, n) .+ 1          # column, row in uncompressed data
    cidx, row = divrem(row-1, step) .+ 1       # codebook index, codebook vector row
    return col, row, cidx
end


Base.setindex!(qa::QuantizedArray, i::Int) =
    @error "setindex! not supported on QuantizedArrays."

Base.setindex!(qa::QuantizedArray, I::Vararg{Int, N}) where {N}=
    @error "setindex! not supported on QuantizedArrays."
