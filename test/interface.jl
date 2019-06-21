@testset "Array interface" begin

    # Test aliases
    @test QuantizedArray(rand(10), k=2, method=:sample) isa QuantizedVector
    @test QuantizedArray(rand(2, 10), k=2, method=:sample) isa QuantizedMatrix

    # Test outer constructor checks
    @test_throws AssertionError QuantizedArray(rand(2, 2, 2), k=2, m=1, method=:sample)
    @test_throws AssertionError QuantizedArray(rand(2, 10), k=-1, m=1, method=:sample)
    @test_throws AssertionError QuantizedArray(rand(2, 10), k=2, m=-1, method=:sample)
    @test_throws AssertionError QuantizedArray(rand(3, 10), k=2, m=2, method=:sample)

    # Test quantize function
    T = Float32
    X = rand(T, 3, 10)
    qa = QuantizedArray(X, k=2, m=3, method=:sample)
    @test QuantizedArrays.quantizer(qa) === qa.quantizer
    @test quantize(QuantizedArrays.quantizer(qa), X) == qa

    # Test other interface parts
    @test size(qa) == qa.quantizer.dims
    @test eltype(qa) == T
    @test QuantizedArrays.nvars(rand(10)) == 1
    @test QuantizedArrays.nvars(rand(2, 10)) == 2

    # Test setindex!
    @test_throws ErrorException qa[1] = one(T)
    @test_throws ErrorException qa[1,1] = one(T)

    # Test getindex

    # vector
    truevector = [0, 1, 0, 2, 0]
    codes = [0x00, 0x01, 0x02]
    vectors = [0 1 2]
    data = [0x00, 0x01, 0x00, 0x02, 0x00]
    cb = CodeBook(codes, vectors)
    aq = ArrayQuantizer(size(truevector), [cb], length(codes), Distances.SqEuclidean())
    qa = QuantizedArray(aq, data)
    @test qa[:] == truevector
    @test all(qa .== truevector)

    # matrix
    truematrix = [0 1 0 2 0 3;
                  0 1 0 2 0 3]
    codes = [0x00, 0x01, 0x02, 0x03]
    vectors = [0 1 2 3;
               0 1 2 3]
    data = [0x00 0x01 0x00 0x02 0x00 0x03]
    cb = CodeBook(codes, vectors)
    aq = ArrayQuantizer(size(truematrix), [cb], length(codes), Distances.SqEuclidean())
    qa = QuantizedArray(aq, data)
    @test qa[:] == truematrix[:]
    @test qa[:,:] == truematrix
    @test all(qa .== truematrix)
end
