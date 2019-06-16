struct ArrayQuantizer{U,D,T,N}
    dims::NTuple{N, Int}                # original array size
    codebooks::Vector{CodeBook{U,D,T}}  # codebooks
    k::Int                              # number of codes/quantizer
end


# show methods
Base.show(io::IO, aq::ArrayQuantizer{U,D,T,N}) where {U,D,T,N} = begin
    m = length(codebooks(aq))
    qstr = ifelse(m==1, "quantizer", "quantizers")
    print(io, "ArrayQuantizer{$U,$D,$T,$N}, $m $qstr, $(aq.k) codes")
end


# Access codebooks
codebooks(aq::ArrayQuantizer) = aq.codebooks


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
                         method::Symbol=DEFAULT_METHOD,
                         distance::Distances.PreMetric=DEFAULT_DISTANCE
                        ) where {T}
    # Get quantized (i.e. key) type
    ncols = size(aa, 2)              # number of columns (samples)
    nrows = size(aa, 1)              # number of rows (variables)
    k = min(k, ncols)                # number of codes for a quantizer
    m = min(m, nrows)                # number of quantizers
    U = quantized_eltype(k)          # type of codes
    D = typeof(distance)

    # Calculate codebooks
    rs = floor(Int, nrows/m)  # row step
    cbooks = Vector{CodeBook{U,D,T}}(undef, m)
    @inbounds @simd for i in 1:m
        rr = rs*(i-1)+1 : rs*i  # row range
        cbooks[i] = build_codebook(aa[rr,:], k, U,
                                   method=method,
                                   distance=distance)
    end
    return ArrayQuantizer(size(aa), cbooks, k)
end

function build_quantizer(aa::AbstractVector{T};
                         k::Int=DEFAULT_K,
                         m::Int=DEFAULT_M,
                         method::Symbol=DEFAULT_METHOD,
                         distance::Distances.PreMetric=DEFAULT_DISTANCE
                        ) where {T}
    aq = build_quantizer(aa', k=k, m=m, method=method, distance=distance)
    return ArrayQuantizer(size(aa), codebooks(aq), k)
end


# Obtain quantized data
function quantize_data(aq::ArrayQuantizer{U,D,T,2}, aa::AbstractMatrix{T}) where {U,D,T}
    nrows, ncols = size(aa)
    @assert nrows == aq.dims[1] "Quantized matrix needs to have $nrows rows"
    cbooks = codebooks(aq)
    m = length(cbooks)
    rs = floor(Int, nrows/m)  # row step
    qa = Matrix{U}(undef, m, ncols)
    @inbounds @simd for i in 1:m
        rr = rs*(i-1)+1 : rs*i  # row range
        qa[i,:] = encode(cbooks[i], aa[rr,:])
    end
    return qa
end

function quantize_data(aq::ArrayQuantizer{U,D,T,1}, aa::AbstractVector{T}) where {U,D,T}
    nrows = aq.dims[1]
    @assert nrows == length(aa) "Quantized vector needs to have $nrows elements"
    aat = aa'  # use transpose as a single row matrix and quantize that
    aqt = ArrayQuantizer(size(aat), codebooks(aq), aq.k)
    qat = quantize_data(aqt, aat)
    return vec(qat)  # return to vector form
end
