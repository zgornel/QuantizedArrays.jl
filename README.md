![Alt text](https://github.com/zgornel/QuantizedArrays.jl/blob/master/docs/src/assets/logo.png)

Array quantization and compression.

[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)
[![Build Status](https://travis-ci.org/zgornel/QuantizedArrays.jl.svg?branch=master)](https://travis-ci.org/zgornel/QuantizedArrays.jl)
[![Coverage Status](https://coveralls.io/repos/github/zgornel/QuantizedArrays.jl/badge.svg?branch=master)](https://coveralls.io/github/zgornel/QuantizedArrays.jl?branch=master)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://zgornel.github.io/QuantizedArrays.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://zgornel.github.io/QuantizedArrays.jl/dev)


## Installation

```julia
Pkg.add("QuantizedArrays")
```


## Examples
 - Vector quantization
```julia
using QuantizedArrays
v = collect(1:10);

# Use 3 quantization levels (always 1 codebook)
qv = QuantizedArray(v, k=3)
# 10-element QuantizedArray{QuantizedArrays.OrthogonalQuantization,UInt8,Distances.SqEuclidean,Int64,1}:
#   5
#   5
#   5
#   5
#   5
#   6
#   6
#   6
#  10
#  10
```

 - matrix quantization
```julia
using QuantizedArrays
m = reshape(collect(1:60), (6,10))
# 6×10 array{int64,2}:
#  1   7  13  19  25  31  37  43  49  55
#  2   8  14  20  26  32  38  44  50  56
#  3   9  15  21  27  33  39  45  51  57
#  4  10  16  22  28  34  40  46  52  58
#  5  11  17  23  29  35  41  47  53  59
#  6  12  18  24  30  36  42  48  54  60

# use 5 quantization levels (vectors) / codebook and 2 codebooks
qm = QuantizedArray(m, k = 5, m = 2)
# 6×10 QuantizedArray{QuantizedArrays.OrthogonalQuantization,UInt8,Distances.SqEuclidean,Int64,2}:
#   7   7  13  13  31  31  43  43  43  55
#   8   8  14  14  32  32  44  44  44  56
#   9   9  15  15  33  33  45  45  45  57
#  10  10  10  28  28  34  40  46  46  46
#  11  11  11  29  29  35  41  47  47  47
#  12  12  12  30  30  36  42  48  48  48
```

## Features
To keep track with the latest features, please consult [NEWS.md](https://github.com/zgornel/QuantizedArrays.jl/blob/master/NEWS.md) and the [documentation](https://zgornel.github.io/QuantizedArrays.jl/dev).


## License

The code has an MIT license and therefore it is free.


## Reporting Bugs

This is work in progress and bugs may still be present...¯\\_(ツ)_/¯ Do not worry, just [open an issue](https://github.com/zgornel/QuantizedArrays.jl/issues/new) to report a bug or request a feature.


## References

 - [Rayuela.jl](https://github.com/una-dinosauria/Rayuela.jl) provides an extensive list of quantization methods implementations as well as good references pertinent to the state-of-the-art.

 - [Quantization in signal processing (Wikipedia)](https://en.wikipedia.org/wiki/Quantization_(signal_processing))

 - [Vector quantization (Wikipedia)](https://en.wikipedia.org/wiki/Vector_quantization)
