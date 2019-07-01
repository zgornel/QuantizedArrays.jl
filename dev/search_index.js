var documenterSearchIndex = {"docs":
[{"location":"#","page":"Introduction","title":"Introduction","text":"CurrentModule=QuantizedArrays","category":"page"},{"location":"#Introduction-1","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"#","page":"Introduction","title":"Introduction","text":"QuantizedArrays is a package for array quantization and compression. The basic principle is that arrays can be represented through the concatenation of shorter vectors obtained by either sampling or clustering of subspaces in the original data. This effectively compresses the original array with a loss in quality.","category":"page"},{"location":"#Implemented-features-1","page":"Introduction","title":"Implemented features","text":"","category":"section"},{"location":"#","page":"Introduction","title":"Introduction","text":"basic AbstractArray interface for vectors, matrices\nquantization types implemented:\nrandom sampling\nk-mean clustering (PQ)\n'cartesian' k-means clustering (OPQ)\nresidual quantization (RVQ)","category":"page"},{"location":"#Installation-1","page":"Introduction","title":"Installation","text":"","category":"section"},{"location":"#","page":"Introduction","title":"Introduction","text":"Installation can be performed from either inside or outside Julia.","category":"page"},{"location":"#Git-cloning-1","page":"Introduction","title":"Git cloning","text":"","category":"section"},{"location":"#","page":"Introduction","title":"Introduction","text":"$ git clone https://github.com/zgornel/QuantizedArrays.jl","category":"page"},{"location":"#Julia-REPL-1","page":"Introduction","title":"Julia REPL","text":"","category":"section"},{"location":"#","page":"Introduction","title":"Introduction","text":"The package can be installed from inside Julia with:","category":"page"},{"location":"#","page":"Introduction","title":"Introduction","text":"using Pkg\nPkg.add(\"QuantizedArrays\")","category":"page"},{"location":"#","page":"Introduction","title":"Introduction","text":"or","category":"page"},{"location":"#","page":"Introduction","title":"Introduction","text":"Pkg.add(PackageSpec(url=\"https://github.com/zgornel/QuantizedArrays.jl\", rev=\"master\"))","category":"page"},{"location":"#","page":"Introduction","title":"Introduction","text":"for the latest master branch.","category":"page"},{"location":"examples/#Usage-examples-1","page":"Usage examples","title":"Usage examples","text":"","category":"section"},{"location":"examples/#","page":"Usage examples","title":"Usage examples","text":"Arrays can be quantized through either the quantize function or using the QuantizedArray constructor directly.","category":"page"},{"location":"examples/#","page":"Usage examples","title":"Usage examples","text":"If the quantization is done through sampling, the quantized arrays may differ with each run as the vector prototypes are randomly sampled.","category":"page"},{"location":"examples/#","page":"Usage examples","title":"Usage examples","text":"using QuantizedArrays\n\nv = collect(1:10)\nqv = quantize(v, k=2)  # quantize by sampling 2 prototypes i.e. values\nqv = QuantizedArray(v, k=2)\n\nm = reshape(collect(1:60), (6,10))\nqm = quantize(m, k=5, m=2)  # 5 prototypes, 2 codebooks\nqm = QuantizedArray(m, k=5, m=2)","category":"page"},{"location":"examples/#","page":"Usage examples","title":"Usage examples","text":"The compressed representation of the input arrays is stored in the data field and the quantizer in the quantizer field","category":"page"},{"location":"examples/#","page":"Usage examples","title":"Usage examples","text":"qm.data\nqm.quantizer","category":"page"},{"location":"examples/#","page":"Usage examples","title":"Usage examples","text":"A new array can be quantized using the quantizers","category":"page"},{"location":"examples/#","page":"Usage examples","title":"Usage examples","text":"quantize(qv.quantizer, rand(1:10, 5))\nquantize(qm.quantizer, rand(1:60, 6, 2))","category":"page"},{"location":"examples/#","page":"Usage examples","title":"Usage examples","text":"The :pq (k-means), :opq ('cartesian' k-means) and :rvq ('residual') quantization methods work for arrays with AbstractFloat elements only and return the same result each run.","category":"page"},{"location":"examples/#","page":"Usage examples","title":"Usage examples","text":"quantize(Float32.(v), k=2, method=:pq)\nquantize(Float32.(m), k=5, m=2, method=:pq)","category":"page"},{"location":"examples/#","page":"Usage examples","title":"Usage examples","text":"Indexing can be performed as in regular arrays","category":"page"},{"location":"examples/#","page":"Usage examples","title":"Usage examples","text":"qm[1,1]\nqm[2,:]\nqm[1:1,:]","category":"page"},{"location":"examples/#","page":"Usage examples","title":"Usage examples","text":"however changing values is not supported","category":"page"},{"location":"examples/#","page":"Usage examples","title":"Usage examples","text":"qm[1,1] = 0","category":"page"},{"location":"examples/#","page":"Usage examples","title":"Usage examples","text":"The element type of the compressed array representation is dynamically calculated to reflect the number of vector prototypes employed","category":"page"},{"location":"examples/#","page":"Usage examples","title":"Usage examples","text":"m = rand(1, 10_000);\nquantize(m, k=256)  # 8 bits are used (UInt8)\nquantize(m, k=257)  # 16 bits are used (UInt16)","category":"page"},{"location":"api/#","page":"API Reference","title":"API Reference","text":"","category":"page"},{"location":"api/#","page":"API Reference","title":"API Reference","text":"Modules = [QuantizedArrays]","category":"page"},{"location":"api/#QuantizedArrays.QuantizedArrays","page":"API Reference","title":"QuantizedArrays.QuantizedArrays","text":"QuantizedArrays.jl - array quantization and compression.\n\n\n\n\n\n","category":"module"},{"location":"api/#QuantizedArrays.ArrayQuantizer","page":"API Reference","title":"QuantizedArrays.ArrayQuantizer","text":"ArrayQuantizer{Q,U,D,T,N}\n\nThe array quantizer object. It transforms an array with elements of type T into a 'quantized' version with elemetns of type U.\n\nFields\n\nquantization::Q type of quantization employed i.e. additive, orthogonal\ndims::NTuple{N,Int} the original array dimensionality\ncodebooks::Vector{CodeBook{U,T}} the codebooks\nk::Int the number of vector prototypes in each codebooks\ndistance::D the distance employed\n\n\n\n\n\n","category":"type"},{"location":"api/#QuantizedArrays.CodeBook","page":"API Reference","title":"QuantizedArrays.CodeBook","text":"CodeBook{U<:Unsigned,T}\n\nThe codebook structure. It holds the codes corresponding to the vector prototypes and the mapping bethween the codes and prototypes.\n\nFields\n\ncodes::Vector{U} the codes\nvectors::Matrix{T} the prototypes\ncodemap::Dict{U,Int} mapping from code to column in vectors\n\n\n\n\n\n","category":"type"},{"location":"api/#QuantizedArrays.QuantizedArray","page":"API Reference","title":"QuantizedArrays.QuantizedArray","text":"QuantizedArray{Q<:AbstractQuantization,U<:Unsigned,D<:PreMetric,T,N} <: AbstractArray{T,N}\n\nA quantized array. It represents a 'quantized' representation of an original array, equivalent to a lossy compressed version.\n\nFields\n\nquantizer::ArrayQuantizer{Q,U,D,T,N} the quantized of the array\ndata::Matrix{U} the actual compressed representation of the quantized array\n\n\n\n\n\n","category":"type"},{"location":"api/#QuantizedArrays.build_codebooks-Union{Tuple{T}, Tuple{U}, Tuple{AbstractArray{T,2},Int64,Int64,Type{U}}} where T where U","page":"API Reference","title":"QuantizedArrays.build_codebooks","text":"build_codebooks(X, k, m, U [;method=DEFAULT_METHOD, distance=DEFAULT_DISTANCE, kwargs])\n\nGenerates m codebooks of k prototypes each for the input matrix X using the algorithm specified my method and distance distance. Specific codebook aglorithm keyword arguments can be specified as well.\n\nArguments\n\nX::AbstractMatrix{T} input matrix of type T\nk::Int number of prototypes/codebook\nm::Int number of codebooks\nU::Type{<:Unsigned} type for codebook codes\n\nKeyword arguments\n\nmethod::Symbol the algorithm to be employed for codebook\n\ngeneration; possible values are :sample (default), :pq for classical k-means clustering codebooks and :opq for 'cartesian' k-means clustering codebooks\n\ndistance::PreMetric the distance to be used in the\n\ncodebook generation methods and data encoding\n\nkwargs... other codebook generation algorithm specific\n\nkeyword arguments such as maxiter::Int.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantizedArrays.build_quantizer-Union{Tuple{AbstractArray}, Tuple{N}, Tuple{T}} where N where T","page":"API Reference","title":"QuantizedArrays.build_quantizer","text":"build_quantizer(aa [;kwargs])\n\nBuilds an array quantizer using the input array aa.\n\nKeyword arguments\n\nk::Int the number of vector prototypes in each codebook\nm::Int the number of codebooks\nmethod::Symbol the algorithm to be employed for codebook\n\ngeneration; possible values are :sample (default), :pq for classical k-means clustering codebooks and :opq for 'cartesian' k-means clustering codebooks\n\ndistance::PreMetric the distance to be used in the\n\ncodebook generation methods and data encoding\n\nOther codebook generation algorithm specific keyword arguments such as maxiter::Int can be specified as well.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantizedArrays.quantize-Union{Tuple{AbstractArray{T,N}}, Tuple{N}, Tuple{T}} where N where T","page":"API Reference","title":"QuantizedArrays.quantize","text":"quantize(aa; [;kwargs])\n\nQuantize an array aa.\n\nArguments\n\naa::AbstractArray input array to be quantized\n\nKeyword arguments\n\nk::Int the number of vector prototypes in each codebook\nm::Int the number of codebooks\nmethod::Symbol the algorithm to be employed for codebook\n\ngeneration; possible values are :sample (default), :pq for classical k-means clustering codebooks and :opq for 'cartesian' k-means clustering codebooks\n\ndistance::PreMetric the distance to be used in the\n\ncodebook generation methods and data encoding\n\nOther codebook generation algorithm specific keyword arguments such as maxiter::Int can be specified as well.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantizedArrays.quantize-Union{Tuple{N}, Tuple{T}, Tuple{D}, Tuple{U}, Tuple{Q}, Tuple{ArrayQuantizer{Q,U,D,T,N},AbstractArray{T,N}}} where N where T where D where U where Q","page":"API Reference","title":"QuantizedArrays.quantize","text":"quantize(aq, aa)\n\nQuantize an array aa using an array quantizer aq.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantizedArrays.AbstractQuantization","page":"API Reference","title":"QuantizedArrays.AbstractQuantization","text":"Abstract quantization type. Its subtypes are two singletons\n\ncorresponding to additive/residual and orthogonal quantization.\n\n\n\n\n\n","category":"type"},{"location":"api/#QuantizedArrays.codebooks-Tuple{ArrayQuantizer}","page":"API Reference","title":"QuantizedArrays.codebooks","text":"codebooks(aq)\n\nAccess the codebooks field of the array quantizer aq.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantizedArrays.decode-Union{Tuple{T}, Tuple{U}, Tuple{CodeBook{U,T},Array{U,1}}} where T where U","page":"API Reference","title":"QuantizedArrays.decode","text":"decode(codebook, codes)\n\nReturns a partial array reconstruction for using the codes and codebook.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantizedArrays.encode-Union{Tuple{T}, Tuple{U}, Tuple{CodeBook{U,T},AbstractArray{T,2}}} where T where U","page":"API Reference","title":"QuantizedArrays.encode","text":"encode(codebook, aa [;distance=DEFAULT_DISTANCE])\n\nEncodes the input array aa using distance to calculate the closest vector prototype from the codebook.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantizedArrays.nvars-Tuple{AbstractArray{T,1} where T}","page":"API Reference","title":"QuantizedArrays.nvars","text":"nvars(aa)\n\nReturns the number of variables of an array aa. For vectors, the value is always 1, for matrices it is the number of rows in the matrix.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantizedArrays.opq_codebooks-Union{Tuple{T}, Tuple{AbstractArray{T,2},Int64,Int64}} where T","page":"API Reference","title":"QuantizedArrays.opq_codebooks","text":"opq_codebooks(X, k, m [;distance=DEFAULT_DISTANCE, maxiter=DEFAULT_PQ_MAXITER])\n\nBuild m codebooks for input matrix X using k sub-space centers obtained using 'cartesian' k-means clustering, the distance distance and maxiter iterations.\n\nReferences:\n\nGe et al. 2014\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantizedArrays.pq_codebooks-Union{Tuple{T}, Tuple{AbstractArray{T,2},Int64,Int64}} where T","page":"API Reference","title":"QuantizedArrays.pq_codebooks","text":"pq_codebooks(X, k, m [;distance=DEFAULT_DISTANCE, maxiter=DEFAULT_PQ_MAXITER])\n\nBuild m codebooks for input matrix X using k sub-space centers obtained using k-means clustering, the distance distance and maxiter iterations.\n\nReferences:\n\nJègou et al. 2011\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantizedArrays.quantization_type-Tuple{Symbol}","page":"API Reference","title":"QuantizedArrays.quantization_type","text":"quantization_type(method)\n\nReturns a quantization type based on the quantization method to be employed.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantizedArrays.quantize_data-Union{Tuple{T}, Tuple{D}, Tuple{U}, Tuple{ArrayQuantizer{OrthogonalQuantization,U,D,T,2},AbstractArray{T,2}}} where T where D where U","page":"API Reference","title":"QuantizedArrays.quantize_data","text":"quantize_data(aq, aa)\n\nReturns a quantized version of the array aa using the array quantizer aq.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantizedArrays.quantized_eltype-Tuple{Any}","page":"API Reference","title":"QuantizedArrays.quantized_eltype","text":"quantized_eltype(k)\n\nDeterminies the minimum Unsigned type that can hold the value k.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantizedArrays.quantizer-Tuple{QuantizedArray}","page":"API Reference","title":"QuantizedArrays.quantizer","text":"quantizer(qa)\n\nAccess the quantizer field of a quantized array qa.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantizedArrays.rowrange-Tuple{Int64,Int64,Int64}","page":"API Reference","title":"QuantizedArrays.rowrange","text":"rowrange(n, m, i)\n\nUtility function that returns a range based on an iteration index i, the number of elements n and number of ranges m.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantizedArrays.rvq_codebooks-Union{Tuple{T}, Tuple{AbstractArray{T,2},Int64,Int64}} where T","page":"API Reference","title":"QuantizedArrays.rvq_codebooks","text":"rvq_codebooks(X, k, m [;distance=DEFAULT_DISTANCE, maxiter=DEFAULT_RVQ_MAXITER])\n\nBuild m codebooks for input matrix X using k layers of residual vector cluster centers, obtained using k-means clustering, the distance distance and maxiter iterations.\n\nReferences:\n\nChen et al. 2010\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantizedArrays.sampling_codebooks-Union{Tuple{T}, Tuple{AbstractArray{T,2},Int64,Int64}} where T","page":"API Reference","title":"QuantizedArrays.sampling_codebooks","text":"sampling_codebooks(X, k, m)\n\nBuild m codebooks for input matrix X by sampling k prototype vectors.\n\n\n\n\n\n","category":"method"}]
}
