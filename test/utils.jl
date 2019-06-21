@testset "Utils" begin
    # rowrange
    n = 23
    v = rand(n)
    for m in 1:9
        step = floor(Int, n/m)
        for i in 1:m
            @test QuantizedArrays.rowrange(n, m, i) ==
                step*(i-1)+1 : step*i
        end
    end

    # quantized_eltype
    for T in subtypes(Unsigned)
        k = rand(typemin(T):typemax(T))
        @test QuantizedArrays.quantized_eltype(k+1) == T
    end
end
