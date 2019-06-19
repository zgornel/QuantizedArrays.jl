# Utility function that returns a row ranged based on an iteration index
@inline rowrange(nrows::Int, m::Int, i::Int) = begin
    rs = floor(Int, nrows/m)  # row step
    rr = rs*(i-1)+1 : rs*i    # row range
    return rr
end

