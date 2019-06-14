struct ArrayQuantizer{U,T,N}
    dims::NTuple{N, Int}                    # original array size
    codebook::Vector{Dict{U, Vector{T}}}    # codes
    k::Int                                  # number of codes/quantizer
end


codebook(aq::ArrayQuantizer{U,T,1}) where {U,T} = aq.codebook[1]

codebook(aq::ArrayQuantizer{U,T,2}) where {U,T} = aq.codebook


# show methods
Base.show(io::IO, aq::ArrayQuantizer{U,T,N}) where {U,T,N} = begin
    m = length(aq.codebook)
    qstr = ifelse(m==1, "quantizer", "quantizers")
    print(io, "ArrayQuantizer{$U,$T,$N}, $m $qstr, $(aq.k) codes")
end


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


_build_sampling_quantizer(aa::AbstractMatrix{T}, k::Int, m::Int) where {T}= begin
    # Get quantized (i.e. key) type
    n = size(aa, 2)                  # number of columns (samples)
    D = size(aa, 1)                  # number of rows (variables)
    k = min(k, n)                    # number of codes for a quantizer
    m = min(m, size(aa,1))           # number of quantizers
    U = quantized_eltype(k)          # type of codes

    # Calculate codebook
    cs = floor(Int, n/k)  # column step
    rs = floor(Int, D/m)  # row step
    cbook = Vector{Dict{U, Vector{T}}}(undef, m)
    @inbounds for j in 1:m
        codes = Dict{U, Vector{T}}()
        rr = rs*(j-1)+1 : rs*j  # row range
        for i in 1:k
            cr = cs*(i-1)+1 : cs*i # column range
            push!(codes, U(i-1)=>aa[rr, rand(cr)])  # randomly sample
        end
        cbook[j] = codes
    end
    return ArrayQuantizer(size(aa), cbook, k)
end

_build_sampling_quantizer(aa::AbstractVector{T}, k::Int, m::Int) where {T}= begin
    q = _build_sampling_quantizer(aa', k, m)
    return ArrayQuantizer(size(aa), codebook(q), k)
end


# Quantization methods
function quantize(aq::ArrayQuantizer{U,T,2}, aa::AbstractMatrix{T}) where {U,T}
    @assert aq.dims == size(aa) "Quantized array needs to have dims=$(aq.dims)."
    cbook = codebook(aq)
    D, n = aq.dims
    m = length(cbook)
    rs = floor(Int, D/m)  # row step
    qa = Matrix{U}(undef, m, n)
    @inbounds for i in 1:m
        rr = rs*(i-1)+1 : rs*i  # row range
        for j in 1:n
            min_dist = Inf
            min_key = zero(U)
            for (k, v) in cbook[i]
                dist = sum((aa[rr,j] .- v).^2)
                dist < min_dist && (min_dist = dist; min_key=k)
            end
            qa[i,j] = min_key
        end
    end
    return qa
end

function quantize(aq::ArrayQuantizer{U,T,1}, aa::AbstractVector{T}) where {U,T}
    @assert aq.dims == size(aa) "Quantized array needs to have dims=$(aq.dims)."
    aat = aa'
    aqt = ArrayQuantizer(size(aat), aq.codebook, aq.k)
    qat = quantize(aqt, aat)
    return vec(qat)
end
