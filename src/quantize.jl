#TODO(Corneliu)

struct ArrayQuantizer{I,T,N}
    dims::NTuple{N, Int}                    # original array size
    codebook::Vector{Dict{I, Vector{T}}}    # codes
    k::Int                                  # number of codes/quantizer
end

codebook(aq::ArrayQuantizer{I,T,1}) where {I,T} = aq.codebook[1]
codebook(aq::ArrayQuantizer{I,T,2}) where {I,T} = aq.codebook


# Quantizer building methods
function build_quantizer(aa::AbstractArray{T,N};
                         k::Int=DEFAULT_K,
                         m::Int=DEFAULT_M,
                         method::Symbol=DEFAULT_METHOD
                        ) where {T,N}
    if method == :sampling
        quantizer = _build_sampling_quantizer(aa, k, m)
    else
        @error "Unknown quantization method '$method'"
    end

    return quantizer
end

# show methods
Base.show(io::IO, aq::ArrayQuantizer{I,T,N}) where {I,T,N} = begin
    m = length(aq.codebook)
    qstr = ifelse(m==1, "quantizer", "quantizers")
    print(io, "ArrayQuantizer{$I,$T,$N}, $m $qstr, $(aq.k) codes")
end


#TODO(Corneliu) implement this
_calculate_quantized_type(k::Int) = begin
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


#TODO Reduce to matrix case i.e. ArrayQuantizer(size(aa), codebook(ArrayQuantizer(aa'),k,m),k)
_build_sampling_quantizer(aa::AbstractVector{T}, k::Int, m::Int) where {T}= begin
    # Get quantized (i.e. key) type
    n = length(aa)
    k = min(k, n)
    I = _calculate_quantized_type(k)

    # Calculate codebook
    d = floor(Int, n/k)
    codes = Dict{I, Vector{T}}()
    @inbounds for i in 1:k
        push!(codes, I(i-1)=>rand(aa[d*(i-1)+1:d*i], 1))  # randomly sample
    end
    #if k*d < n
    #    push!(codes, I(k-1)=>rand(aa[k*d+1:n], 1))  # last elements
    #end
    return ArrayQuantizer(size(aa), [codes], k)
end


_build_sampling_quantizer(aa::AbstractMatrix{T}, k::Int, m::Int) where {T}= begin
    # Get quantized (i.e. key) type
    n = size(aa, 2)                  # number of columns (samples)
    D = size(aa, 1)                  # number of rows (variables)
    k = min(k, n)                    # number of codes for aquantizer
    m = min(m, size(aa,1))           # number of quantizers
    I = _calculate_quantized_type(k) # type of codes

    # Calculate codebook
    cs = floor(Int, n/k)  # column step
    rs = floor(Int, D/m)  # row step
    cbook = Vector{Dict{I, Vector{T}}}(undef, m)
    @inbounds for j in 1:m
        codes = Dict{I, Vector{T}}()
        rr = rs*(j-1)+1 : rs*j  # row range
        for i in 1:k
            cr = cs*(i-1)+1 : cs*i # column range
            push!(codes, I(i-1)=>aa[rr, rand(cr)])  # randomly sample
        end
        cbook[j] = codes
    end
    return ArrayQuantizer(size(aa), cbook, k)
end


# Quantization methods
# TODO Reduce to matrix form, if possible
function quantize(aq::ArrayQuantizer{I,T,1}, aa::AbstractVector{T}) where {I,T}
    @assert aq.dims == size(aa) "ArrayQuantizer dimension mismatch"
    qa = Vector{I}(undef, aq.dims...)
    cbook = codebook(aq)

    @inbounds for i in 1:length(aa)
        min_dist = Inf
        min_key = zero(I)
        for (k, v) in cbook
            dist = (aa[i] - v[1])^2
            dist < min_dist && (min_dist = dist; min_key=k)
        end
        qa[i] = min_key
    end
    return qa
end

function quantize(aq::ArrayQuantizer{I,T,2}, aa::AbstractMatrix{T}) where {I,T}
    @assert aq.dims == size(aa) "ArrayQuantizer dimension mismatch"
    cbook = codebook(aq)
    D, n = aq.dims
    m = length(cbook)
    rs = floor(Int, D/m)  # row step
    qa = Matrix{I}(undef, m, n)
    @inbounds for i in 1:m
        rr = rs*(i-1)+1 : rs*i  # row range
        for j in 1:n
            min_dist = Inf
            min_key = zero(I)
            for (k, v) in cbook[i]
                dist = sum((aa[rr,j] .- v).^2)
                dist < min_dist && (min_dist = dist; min_key=k)
            end
            qa[i,j] = min_key
        end
    end
    return qa
end
