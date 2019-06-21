# Utility function that returns a row ranged based on an iteration index
@inline rowrange(nrows::Int, m::Int, i::Int) = begin
    rs = floor(Int, nrows/m)  # row step
    rr = rs*(i-1)+1 : rs*i    # row range
    return rr
end


# Determine code type
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
