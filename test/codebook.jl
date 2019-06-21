@testset "Codebook generation" begin

# Generate some random data
n = 100
dim = 2
U = UInt8
Xfloat = rand(dim, n)
codes = U.(collect(1:n))

# Test constructors
cb = CodeBook(codes, Xfloat)
@test cb isa CodeBook{U, eltype(Xfloat)}
@test size(cb.vectors) == size(Xfloat)
@test length(cb.codes) == length(codes)
for (code, idx) in cb.codemap  # tests codemap creation, getindex
    @test cb[code] == Xfloat[:, cb.codemap[code]] == Xfloat[:, idx]
end
@test_throws DimensionMismatch CodeBook(codes, rand(2,n-3))
@test_throws ErrorException CodeBook{UInt8}(rand(2,257))

for T in [UInt8, UInt16, UInt32, UInt64, Int8, Int16, Int32, Int64,
          Float16, Float32, Float64]
    for U in [UInt8, UInt16]  # not really necessary for UInt32, UInt64
        X = rand(T, dim, Int.(typemax(U)))
        @test CodeBook(X) isa CodeBook{U,T}
        if T<:Unsigned && QuantizedArrays.TYPE_TO_BITS[T] <= QuantizedArrays.TYPE_TO_BITS[U]
            @test CodeBook{U,U}(X) isa CodeBook{U,U}
        elseif T<:Unsigned && QuantizedArrays.TYPE_TO_BITS[T] >= QuantizedArrays.TYPE_TO_BITS[U]
            @test_throws InexactError CodeBook{U,U}(X)
        elseif T<:AbstractFloat
            @test CodeBook{U, Float64}(X) isa CodeBook{U, Float64}
        end
    end
end

@test QuantizedArrays.encode(cb, Xfloat) == cb.codes

#TODO(Corneliu): Functional testing of codebook
#                generation algorithms i.e. `build_codebooks`

end
