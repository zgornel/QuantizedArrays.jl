@testset "Array interface" begin

    # Test aliases
    @test QuantizedArray(rand(10), k=2, method=:sample) isa QuantizedVector
    @test QuantizedArray(rand(2, 10), k=2, method=:sample) isa QuantizedMatrix

    # Test outer constructor checks
    @test_throws AssertionError QuantizedArray(rand(2, 2, 2), k=2, m=1, method=:sample)
    @test_throws AssertionError QuantizedArray(rand(2, 10), k=-1, m=1, method=:sample)
    @test_throws AssertionError QuantizedArray(rand(2, 10), k=2, m=-1, method=:sample)
    @test_throws AssertionError QuantizedArray(rand(3, 10), k=2, m=2, method=:sample)
    for method in [:pq, :opq, :rvq]
        @test_throws AssertionError QuantizedArray(rand([1,2,3], 3, 10), k=2, m=2, method=method)
    end

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

    # Test orthogonal reconstruction: getindex

    # vector
    truevector = [0, 1, 0, 2, 0]
    codes = [0x00, 0x01, 0x02]
    vectors = [0 1 2]
    data = [0x00 0x01 0x00 0x02 0x00]
    cb = CodeBook(codes, vectors)
    q = QuantizedArrays.OrthogonalQuantization()
    rot = diagm(0=>ones(eltype(truevector), length(truevector)))
    aq = ArrayQuantizer(q, size(truevector), [cb], length(codes), Distances.SqEuclidean(), rot)
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
    q = QuantizedArrays.OrthogonalQuantization()
    rot = diagm(0=>ones(eltype(truematrix), size(truematrix, 1)))
    aq = ArrayQuantizer(q, size(truematrix), [cb], length(codes), Distances.SqEuclidean(), rot)
    qa = QuantizedArray(aq, data)
    @test qa[:] == truematrix[:]
    @test qa[:,:] == truematrix
    @test all(qa .== truematrix)

    # Test aditive reconstruction: getindex

    # vector
    truevector = [0, 1, 0, 2, 0] .+ [-1, -1, 1, 1, 2]
    codes = [0x00, 0x01, 0x02]
    vectors1 = [0 1 2]
    vectors2 = [-1 1 2]
    data = [0x00 0x01 0x00 0x02 0x00;
            0x00 0x00 0x01 0x01 0x02]
    cbs = [CodeBook(codes, vectors1), CodeBook(codes, vectors2)]
    q = QuantizedArrays.AdditiveQuantization()
    rot = diagm(0=>ones(eltype(truevector), length(truevector)))
    aq = ArrayQuantizer(q, size(truevector), cbs, length(codes), Distances.SqEuclidean(), rot)
    qa = QuantizedArray(aq, data)
    @test qa[:] == truevector
    @test all(qa .== truevector)

    # matrix
    truematrix = [0 1 0 2 0 3;
                  0 1 0 2 0 3] .+
                 [ 2  3  1  4  4  1;
                  -2 -3 -1 -4 -4 -1]
    codes = [0x00, 0x01, 0x02, 0x03]
    vectors1 = [0 1 2 3;
                0 1 2 3]
    vectors2 = [1  2  3  4
                -1 -2 -3 -4]
    data = [0x00 0x01 0x00 0x02 0x00 0x03;
            0x01 0x02 0x00 0x03 0x03 0x00]
    cbs = [CodeBook(codes, vectors1), CodeBook(codes, vectors2)]
    q = QuantizedArrays.AdditiveQuantization()
    rot = diagm(0=>ones(eltype(truematrix), size(truematrix, 1)))
    aq = ArrayQuantizer(q, size(truematrix), cbs, length(codes), Distances.SqEuclidean(), rot)
    qa = QuantizedArray(aq, data)
    @test qa[:] == truematrix[:]
    @test qa[:,:] == truematrix
    @test all(qa .== truematrix)

    # Setindex tests
    @test_throws ErrorException qa[1] = 1
    @test_throws ErrorException qa[1, 1] = 1

end
