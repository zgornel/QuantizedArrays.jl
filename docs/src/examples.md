# Usage examples

Arrays can be quantized through either the `quantize` function or using the `QuantizedArray` constructor directly.

If the quantization is done through sampling, the quantized arrays may differ with each run as the vector prototypes are randomly sampled.
```@repl index
using QuantizedArrays

v = collect(1:10)
qv = quantize(v, k=2)  # quantize by sampling 2 prototypes i.e. values
qv = QuantizedArray(v, k=2)

m = reshape(collect(1:60), (6,10))
qm = quantize(m, k=5, m=2)  # 5 prototypes, 2 codebooks
qm = QuantizedArray(m, k=5, m=2)
```

The compressed representation of the input arrays is stored in the `data` field and the quantizer in the `quantizer` field
```@repl index
qm.data
qm.quantizer
```
A new array can be quantized using the quantizers
```@repl index
quantize(qv.quantizer, rand(1:10, 5))
quantize(qm.quantizer, rand(1:60, 6, 2))
```

The `:pq` (k-means), `:opq` ('cartesian' k-means) and `:rvq` ('residual') quantization methods work for arrays with `AbstractFloat` elements only and return the same result each run.
```@repl index
quantize(Float32.(v), k=2, method=:pq)
quantize(Float32.(m), k=5, m=2, method=:pq)
```

Indexing can be performed as in regular arrays
```@repl index
qm[1,1]
qm[2,:]
qm[1:1,:]
```
however changing values is not supported
```@repl index
qm[1,1] = 0
```

The element type of the compressed array representation is dynamically calculated to reflect the number of vector prototypes employed
```@repl index
m = rand(1, 10_000);
quantize(m, k=256)  # 8 bits are used (UInt8)
quantize(m, k=257)  # 16 bits are used (UInt16)
```
