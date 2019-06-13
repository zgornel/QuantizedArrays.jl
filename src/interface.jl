struct QuantizedArray{I<:Unsigned, T, N} <: AbstractArray{T,N}
    data::Array{I, N}                       # compressed data
    quantizer::ArrayQuantizer{I,T,N}        # codebook: vector prototypes
end


# Aliases
const QuantizedVector{I,T} = QuantizedArray{I,T,1}
const QuantizedMatrix{I,T} = QuantizedArray{I,T,2}


# Basic constructor
function QuantizedArray(aa::AbstractArray{T,N};
                        k::Int=DEFAULT_K,
                        m::Int=DEFAULT_M,
                        method::Symbol=DEFAULT_METHOD,
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
# Base.show(io::IO, ::MIME"text/plain", qa::QuantizedVector{I,T}) where {I,T} = begin
#     print(io, "$(length(qa))-element QuantizedArray{$I,$T,1}")
#     print(io, "…")
# end
#
# Base.show(io::IO, ::MIME"text/plain", qa::QuantizedMatrix{I,T}) where {I,T} = begin
#     print(io, "$(qa.quantizer.dims[1])×$(qa.quantizer.dims[2]) QuantizedArray{$I,$T,2}")
#     print(io, "…")
# end


# eltype, size, length
Base.eltype(qa::QuantizedArray{I,T,N}) where {I,T,N} = T

Base.size(qa::QuantizedArray) = qa.quantizer.dims

Base.length(qa::QuantizedArray) = prod(qa.quantizer.dims)


### # Similar methods
### #TODO(Corneliu) Other `similar` emthods ?
### function Base.similar(qa::QuantizedArray{I,T,N}) where {I,T,N}
###     _data = similar(qa.data)
###     _codebook = similar(qa.codebook)
###     return QuantizedArray(qa.D, qa.k, qa.m, qa.dims, _data, _codebook)
### end


# Indexing interface: getindex, setindex!
Base.IndexStyle(::Type{QuantizedArray}) = IndexLinear()

@inline function Base.getindex(qa::QuantizedVector, i::Int)
    @boundscheck checkbounds(qa.data, i)
    #qi = floor(Int, i/qa.k)  # corresponding quantized index
    qkey = getindex(qa.data, i)  # get quantization code
    cb = codebook(qa.quantizer)  # get codebook
    return cb[qkey][1]           # get quantized value
end

@inline function Base.getindex(qa::QuantizedMatrix, i::Int)
    cbook = codebook(quantizer(qa))
    m = length(cbook)
    D, n = size(qa)
    iq, ii = divrem(i-1, Int(D/m))  # quantized index, partial index
    qkey = getindex(qa.data, iq+1)  # get quantization code
    return cbook[iq+1][qkey][ii+1]  # get quantized value
end

@inline function Base.getindex(qa::QuantizedMatrix, i::Int, j::Int)
    cbook = codebook(quantizer(qa))
    m = length(cbook)
    D, n = size(qa)
    rs = floor(Int, D/m)               # row step
    iq, ii = divrem(i-1, rs)           # quantized index, partial index
    qkey = getindex(qa.data, iq+1, j)  # get quantization code
    return cbook[iq+1][qkey][ii+1]     # get quantized value
end
