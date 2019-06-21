@testset "Array quantizer" begin

m = 2
k = 2
for T in [UInt8, UInt16, UInt32, UInt64, Int8, Int16, Int32, Int64,
          Float16, Float32, Float64]
    for mᵢ in 1:2
        X = rand(T[0,1,2,3], m, 10)  # values in 0-3 to keep distances small
        aqm = build_quantizer(X, k=k, m=mᵢ, method=:sample)
        @test aqm isa ArrayQuantizer{UInt8, Distances.SqEuclidean, T, 2}
        @test length(aqm.codebooks) == mᵢ
        qm = QuantizedArrays.quantize_data(aqm, X)
        @test qm isa Matrix{UInt8}
        @test size(qm) == (mᵢ, size(X,2))

        Xᵥ= rand(T[0,1,2,3,4], 10)
        aqv = build_quantizer(Xᵥ, k=k, m=mᵢ, method=:sample)
        @test aqv isa ArrayQuantizer{UInt8, Distances.SqEuclidean, T, 1}
        @test length(aqv.codebooks) == 1
        qv = QuantizedArrays.quantize_data(aqv, Xᵥ)
        @test qv isa Vector{UInt8}
        @test length(qv) == length(Xᵥ)
    end
end


T = Float32
Twrong = Float16
X = rand(T, m, 100)
m = 2
k = 2
aq = build_quantizer(X, k=k, m=m, method=:sample)
@test QuantizedArrays.codebooks(aq) == aq.codebooks

@test_throws AssertionError QuantizedArrays.quantize_data(aq, rand(T, m+1, 100))
@test_throws MethodError QuantizedArrays.quantize_data(aq, rand(Twrong, m, 100))

# The ArrayQuantizer works with any `m`:
X2 = rand(T, m+1, 100)
aq2 = ArrayQuantizer(X2, k=k, m=m, method=:sample)
@test aq2 isa ArrayQuantizer{UInt8, Distances.SqEuclidean, T, 2}
@test length(aq2.codebooks) == m

aq3 = ArrayQuantizer(X, k=k, m=m+100, method=:sample)
@test aq3 isa ArrayQuantizer{UInt8, Distances.SqEuclidean, T, 2}
@test length(aq3.codebooks) == m

end
