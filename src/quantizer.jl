struct ArrayQuantizer{U,D,T,N}
    dims::NTuple{N, Int}              # original array size
    codebooks::Vector{CodeBook{U,T}}  # codebooks
    k::Int                            # number of codes/quantizer
    distance::D
end


# show methods
Base.show(io::IO, aq::ArrayQuantizer{U,D,T,N}) where {U,D,T,N} = begin
    m = length(codebooks(aq))
    qstr = ifelse(m==1, "quantizer", "quantizers")
    cstr = ifelse(aq.k==1, "code", "codes")
    print(io, "ArrayQuantizer{$U,$D,$T,$N}, $m $qstr, $(aq.k) $cstr")
end


# Access codebooks
codebooks(aq::ArrayQuantizer) = aq.codebooks


# Quantizer building methods
function build_quantizer(aa::AbstractMatrix{T};
                         k::Int=DEFAULT_K,
                         m::Int=DEFAULT_M,
                         method::Symbol=DEFAULT_METHOD,
                         distance::Distances.PreMetric=DEFAULT_DISTANCE
                        ) where {T}
    cbooks = build_codebooks(aa, k, m, method=method, distance=distance)
    return ArrayQuantizer(size(aa), cbooks, k, distance)
end

function build_quantizer(aa::AbstractVector{T};
                         k::Int=DEFAULT_K,
                         m::Int=DEFAULT_M,
                         method::Symbol=DEFAULT_METHOD,
                         distance::Distances.PreMetric=DEFAULT_DISTANCE
                        ) where {T}
    aq = build_quantizer(aa', k=k, m=m, method=method, distance=distance)
    return ArrayQuantizer(size(aa), codebooks(aq), k, distance)
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
        qa[i,:] = encode(cbooks[i], aa[rr,:], distance=aq.distance)
    end
    return qa
end

function quantize_data(aq::ArrayQuantizer{U,D,T,1}, aa::AbstractVector{T}) where {U,D,T}
    nrows = aq.dims[1]
    @assert nrows == length(aa) "Quantized vector needs to have $nrows elements"
    aat = aa'  # use transpose as a single row matrix and quantize that
    aqt = ArrayQuantizer(size(aat), codebooks(aq), aq.k, aq.distance)
    qat = quantize_data(aqt, aat)
    return vec(qat)  # return to vector form
end
