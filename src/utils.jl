"""
    rowrange(n, m, i)

Utility function that returns a range based on an iteration index `i`,
the number of elements `n` and number of ranges `m`.
"""
@inline rowrange(n::Int, m::Int, i::Int) = begin
    rs = floor(Int, n/m)  # row step
    rr = rs*(i-1)+1 : rs*i    # row range
    return rr
end


"""
    quantized_eltype(k)

Determinies the minimum `Unsigned` type that can hold the value `k`.
"""
quantized_eltype(k) = begin
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


"""
    quantization_type(method)

Returns a quantization type based on the quantization method to be employed.
"""
quantization_type(method::Symbol) = begin
    ifelse(method == :rvq, AdditiveQuantization(), OrthogonalQuantization())
end
