struct QuantizedArray{U<:Unsigned, T, N} <: AbstractArray{T,N}
    data::Array{U, N}                       # compressed data
    quantizer::ArrayQuantizer{U,T,N}        # codebook: vector prototypes
end


# Aliases
const QuantizedVector{U,T} = QuantizedArray{U,T,1}
const QuantizedMatrix{U,T} = QuantizedArray{U,T,2}


# Basic constructor
function QuantizedArray(aa::AbstractArray{T,N};
                        k::Int=DEFAULT_K,
                        m::Int=DEFAULT_M,
                        method::Symbol=DEFAULT_METHOD
                       ) where {T,N}
    @assert N <=2 "Array quantization is supported only for Vectors and Matrices"
    @assert k >= 1 "`k` has to be larger or equal to 1"
    @assert m >= 1 "`m` has to be larger or equal to 1"
    if N == 2
        nrows = size(aa, 1)
        @assert rem(nrows, m) == 0 "`m` has to divide exactly $nrows rows"
    end
    quantizer = build_quantizer(aa, k=k, m=m, method=method)
    qa = quantize(quantizer, aa)
    return QuantizedArray(qa, quantizer)
end


# Get quantizer function
quantizer(qa::QuantizedArray) = qa.quantizer

# show methods
# Base.show(io::IO, ::MIME"text/plain", qa::QuantizedVector{U,T}) where {U,T} = begin
#     print(io, "$(length(qa))-element QuantizedArray{$U,$T,1}")
#     print(io, "…")
# end
#
# Base.show(io::IO, ::MIME"text/plain", qa::QuantizedMatrix{U,T}) where {U,T} = begin
#     print(io, "$(qa.quantizer.dims[1])×$(qa.quantizer.dims[2]) QuantizedArray{$U,$T,2}")
#     print(io, "…")
# end


# eltype, size, length
Base.eltype(qa::QuantizedArray{U,T,N}) where {U,T,N} = T

Base.size(qa::QuantizedArray) = qa.quantizer.dims

Base.length(qa::QuantizedArray) = prod(qa.quantizer.dims)


### # Similar methods
### #TODO(Corneliu) Other `similar` emthods ?
### function Base.similar(qa::QuantizedArray{U,T,N}) where {U,T,N}
###     _data = similar(qa.data)
###     _codebook = similar(qa.codebook)
###     return QuantizedArray(qa.D, qa.k, qa.m, qa.dims, _data, _codebook)
### end

Base.IndexStyle(::Type{<:QuantizedArray}) = IndexLinear()


# Indexing interface: getindex, setindex!
@inline function Base.getindex(qa::QuantizedVector, i::Int)
    @boundscheck checkbounds(qa.data, i)
    qkey = getindex(qa.data, i)  # get quantization code
    cb = codebook(qa.quantizer)  # get codebook
    return cb[qkey][1]           # get quantized value
end

@inline function Base.getindex(qa::QuantizedMatrix, i::Int)
    cbook = codebook(quantizer(qa))
    m = length(cbook)
    D, n = size(qa)
    iq, ii = _quantized_indices(Int(D/m), I...)  # quantized index, partial index
    qkey = getindex(qa.data, iq)  # get quantization code
    return cbook[iq][qkey][ii]    # get quantized value
end

@inline function Base.getindex(qa::QuantizedMatrix, i::Int, j::Int)
    cbook = codebook(quantizer(qa))
    m = length(cbook)
    D, n = size(qa)
    iq, ii = _quantized_indices(Int(D/m), i, j)  # quantized index, partial index
    qkey = getindex(qa.data, iq, j)  # get quantization code
    return cbook[iq][qkey][ii]          # get quantized value
end

function _quantized_indices(step::Int, I::Vararg{Int,N}) where {N}
    iq, ii = divrem(I[1]-1, step) .+ 1
    return iq, ii
end
