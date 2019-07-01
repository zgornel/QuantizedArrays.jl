"""
    QuantizedArray{Q<:AbstractQuantization,U<:Unsigned,D<:PreMetric,T,N} <: AbstractArray{T,N}

A quantized array. It represents a 'quantized' representation of an original
array, equivalent to a lossy compressed version.

# Fields
  * `quantizer::ArrayQuantizer{Q,U,D,T,N}` the quantized of the array
  * `data::Matrix{U}` the actual compressed representation of the quantized array
"""
struct QuantizedArray{Q,U<:Unsigned,D<:Distances.PreMetric,T,N} <: AbstractArray{T,N}
    quantizer::ArrayQuantizer{Q,U,D,T,N}
    data::Matrix{U}
end


# Aliases
const QuantizedVector{Q,U,D,T} = QuantizedArray{Q,U,D,T,1}
const QuantizedMatrix{Q,U,D,T} = QuantizedArray{Q,U,D,T,2}
const OrthogonalQuantizedArray{U,D,T,N} = QuantizedArray{OrthogonalQuantization,U,D,T,N}
const OrthogonalQuantizedMatrix{U,D,T,N} = QuantizedArray{OrthogonalQuantization,U,D,T,2}
const AdditiveQuantizedArray{U,D,T,N} = QuantizedArray{AdditiveQuantization,U,D,T,N}
const AdditiveQuantizedMatrix{U,D,T,N} = QuantizedArray{AdditiveQuantization,U,D,T,2}


# Main outer constructor
function QuantizedArray(aa::AbstractArray{T,N};
                        k::Int=DEFAULT_K,
                        m::Int=DEFAULT_M,
                        method::Symbol=DEFAULT_METHOD,
                        distance::Distances.PreMetric=DEFAULT_DISTANCE,
                        kwargs...) where {T,N}
    @assert N <=2 "Array quantization is supported only for Vectors and Matrices"
    @assert k >= 1 "`k` has to be larger or equal to 1"
    @assert m >= 1 "`m` has to be larger or equal to 1"
    method != :rvq && @assert rem(nvars(aa), m) == 0 "`m` has to divide exactly $(nvars(aa))"
    aq = build_quantizer(aa, k=k, m=m, method=method, distance=distance; kwargs...)
    data = quantize_data(aq, aa)
    return QuantizedArray(aq, data)
end


"""
    quantize(aq, aa)

Quantize an array `aa` using an array quantizer `aq`.
"""
function quantize(aq::ArrayQuantizer{Q,U,D,T,N}, aa::AbstractArray{T,N}) where {Q,U,D,T,N}
    new_aq = ArrayQuantizer(aq.quantization, size(aa), codebooks(aq), aq.k, aq.distance)
    data = quantize_data(new_aq, aa)
    return QuantizedArray(new_aq, data)
end


"""
    quantize(aa; [;kwargs])
Quantize an array `aa`.

# Arguments
  * `aa::AbstractArray` input array to be quantized

# Keyword arguments
  * `k::Int` the number of vector prototypes in each codebook
  * `m::Int` the number of codebooks
  * `method::Symbol` the algorithm to be employed for codebook
generation; possible values are `:sample` (default), `:pq` for
classical k-means clustering codebooks and `:opq` for 'cartesian'
k-means clustering codebooks
  * `distance::PreMetric` the distance to be used in the
codebook generation methods and data encoding

Other codebook generation algorithm specific keyword arguments
such as `maxiter::Int` can be specified as well.
"""
quantize(aa::AbstractArray{T,N}; kwargs...) where {T,N} = QuantizedArray(aa; kwargs...)


# eltype, size, length
Base.eltype(qa::QuantizedArray{Q,U,D,T,N}) where {Q,U,D,T,N} = T
Base.size(qa::QuantizedArray) = qa.quantizer.dims
Base.IndexStyle(::Type{<:QuantizedArray}) = IndexLinear()


"""
    quantizer(qa)

Access the quantizer field of a quantized array `qa`.
"""
quantizer(qa::QuantizedArray) = qa.quantizer


"""
    nvars(aa)

Returns the number of variables of an array `aa`. For vectors,
the value is always `1`, for matrices it is the number of
rows in the matrix.
"""
nvars(av::AbstractVector) = 1
nvars(am::AbstractMatrix) = size(am, 1)


# Indexing interface: getindex, setindex!
@inline function Base.getindex(qa::OrthogonalQuantizedArray, i::Int)
    cbooks = codebooks(quantizer(qa))
    m = length(cbooks)
    col, row, cidx = _indices(nvars(qa), m, i)
    qkey = getindex(qa.data, cidx, col)
    return cbooks[cidx][qkey][row]
end

@inline function Base.getindex(qa::OrthogonalQuantizedMatrix, i::Int, j::Int)
    cbooks = codebooks(quantizer(qa))
    m = length(cbooks)
    _, row, cidx = _indices(nvars(qa), m, i)
    qkey = getindex(qa.data, cidx, j)
    return cbooks[cidx][qkey][row]
end

@inline function Base.getindex(qa::AdditiveQuantizedArray, i::Int)
    cbooks = codebooks(quantizer(qa))
    m = length(cbooks)
    col, row = divrem(i[1]-1, nvars(qa)) .+ 1
    return sum(cbooks[cidx][getindex(qa.data, cidx, col)][row] for cidx in 1:m)
end

@inline function Base.getindex(qa::AdditiveQuantizedMatrix, i::Int, j::Int)
    cbooks = codebooks(quantizer(qa))
    m = length(cbooks)
    return sum(cbooks[cidx][getindex(qa.data, cidx, j)][i] for cidx in 1:m)
end

function _indices(n::Int, m::Int, I::Vararg{Int,N}) where {N}
    step = Int(n/m)
    col, row = divrem(I[1]-1, n) .+ 1          # column, row in uncompressed data
    cidx, row = divrem(row-1, step) .+ 1       # codebook index, codebook vector row
    return col, row, cidx
end


# setindex! not supported
Base.setindex!(qa::QuantizedArray, i::Int) =
    @error "setindex! not supported on QuantizedArrays"

Base.setindex!(qa::QuantizedArray, I::Vararg{Int, N}) where {N}=
    @error "setindex! not supported on QuantizedArrays"
