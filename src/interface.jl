struct QuantizedArray{U<:Unsigned,D<:Distances.PreMetric,T,N} <: AbstractArray{T,N}
    data::Array{U,N}
    quantizer::ArrayQuantizer{U,D,T,N}
end


# Aliases
const QuantizedVector{U,D,T} = QuantizedArray{U,D,T,1}
const QuantizedMatrix{U,D,T} = QuantizedArray{U,D,T,2}


# Main outer constructor


function QuantizedArray(aa::AbstractArray{T,N};
                        k::Int=DEFAULT_K,
                        m::Int=DEFAULT_M,
                        method::Symbol=DEFAULT_METHOD,
                        distance::Distances.PreMetric=DEFAULT_DISTANCE) where {T,N}
    @assert N <=2 "Array quantization is supported only for Vectors and Matrices"
    @assert k >= 1 "`k` has to be larger or equal to 1"
    @assert m >= 1 "`m` has to be larger or equal to 1"
    @assert rem(nvars(aa), m) == 0 "`m` has to divide exactly $(nvars(aa))"
    aq = build_quantizer(aa, k=k, m=m, method=method, distance=distance)
    data = quantize_data(aq, aa)
    return QuantizedArray(data, aq)
end


# Quantize an array based on an external codebook
function quantize(aq::ArrayQuantizer{U,D,T,N},
                  aa::AbstractArray{T,N}) where {U,D,T,N}
    new_aq = ArrayQuantizer(size(aa), codebooks(aq), aq.k, aq.distance)
    data = quantize_data(new_aq, aa)
    return QuantizedArray(data, new_aq)
end

quantize(aa::AbstractArray{T,N}; kwargs...) where {T,N} = quantize(aa; kwargs...)


# eltype, size, length
Base.eltype(qa::QuantizedArray{U,D,T,N}) where {U,D,T,N} = T
Base.size(qa::QuantizedArray) = qa.quantizer.dims
Base.IndexStyle(::Type{<:QuantizedArray}) = IndexLinear()


# Access the quantizer
quantizer(qa::QuantizedArray) = qa.quantizer


# nvars:
#   - vectors have a single variable
#   - matrices have number of rows variables
nvars(av::AbstractVector) = 1
nvars(am::AbstractMatrix) = size(am, 1)


# Indexing interface: getindex, setindex!
@inline function Base.getindex(qa::QuantizedArray, i::Int)
    cbooks = codebooks(quantizer(qa))
    m = length(cbooks)
    col, row, cidx = _indices(nvars(qa), m, i)   # get various indexing numbers
    qkey = getindex(qa.data, m*(col-1)+cidx)     # get quantization code
    return cbooks[cidx][qkey][row]               # get quantized value
end

@inline function Base.getindex(qa::QuantizedMatrix, i::Int, j::Int)
    cbooks = codebooks(quantizer(qa))
    m = length(cbooks)
    _, row, cidx = _indices(nvars(qa), m, i)   # get various indexing numbers
    qkey = getindex(qa.data, cidx, j)          # get quantization code
    return cbooks[cidx][qkey][row]             # get quantized value
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
