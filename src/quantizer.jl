struct ArrayQuantizer{U,T,N}
    dims::NTuple{N, Int}              # original array size
    codebooks::Vector{CodeBook{U,T}}  # codebooks
    k::Int                            # number of codes/quantizer
end


# Access codebooks
codebooks(aq::ArrayQuantizer{U,T,1}) where {U,T} = aq.codebooks[1]
codebooks(aq::ArrayQuantizer{U,T,2}) where {U,T} = aq.codebooks


# show methods
Base.show(io::IO, aq::ArrayQuantizer{U,T,N}) where {U,T,N} = begin
    m = length(aq.codebooks)
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
function build_quantizer(aa::AbstractMatrix{T};
                         k::Int=DEFAULT_K,
                         m::Int=DEFAULT_M,
                         method::Symbol=DEFAULT_METHOD) where {T}
    # Get quantized (i.e. key) type
    n = size(aa, 2)                  # number of columns (samples)
    D = size(aa, 1)                  # number of rows (variables)
    k = min(k, n)                    # number of codes for a quantizer
    m = min(m, size(aa,1))           # number of quantizers
    U = quantized_eltype(k)          # type of codes

    # Calculate codebooks
    rs = floor(Int, D/m)  # row step
    cbooks = Vector{CodeBook{U,T}}(undef, m)
    @inbounds @simd for i in 1:m
        rr = rs*(i-1)+1 : rs*i  # row range
        cbooks[i] = build_codebook(aa[rr,:], k, U, method=method)
    end
    return ArrayQuantizer(size(aa), cbooks, k)
end

function build_quantizer(aa::AbstractVector{T};
                         k::Int=DEFAULT_K,
                         m::Int=DEFAULT_M,
                         method::Symbol=DEFAULT_METHOD) where {T}
    q = build_quantizer(aa', k=k, m=m, method=method)
    return ArrayQuantizer(size(aa), codebooks(q), k)
end


# Quantization methods
function quantize(aq::ArrayQuantizer{U,T,2},
                  aa::AbstractMatrix{T};
                  distance=Distances.Euclidean()) where {U,T}
    @assert aq.dims == size(aa) "Quantized array needs to have dims=$(aq.dims)."
    cbooks = codebooks(aq)
    D, n = aq.dims
    m = length(cbooks)
    rs = floor(Int, D/m)  # row step
    qa = Matrix{U}(undef, m, n)
    @inbounds @simd for i in 1:m
        rr = rs*(i-1)+1 : rs*i  # row range
        qa[i,:] = infer_codes(cbooks[i], aa[rr,:], distance=distance)
    end
    return qa
end

function quantize(aq::ArrayQuantizer{U,T,1},
                  aa::AbstractVector{T};
                  distance=Distances.Euclidean()) where {U,T}
    @assert aq.dims == size(aa) "Quantized array needs to have dims=$(aq.dims)."
    aat = aa'  # use transpose as a single row matrix and quantize that
    aqt = ArrayQuantizer(size(aat), aq.codebooks, aq.k)
    qat = quantize(aqt, aat, distance=distance)
    return vec(qat)  # return to vector form
end
