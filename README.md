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


## TODO
~~ - basic interface for vectors, matrices (w. sampling quantization) ~~
 - extend basic array interface
 - actual quantization (i.e. cluster sub-vectors w. PQ/OPQ/etc.)
 - tests
 - documetation
 - implement higher order tensor interface


## Examples
```julia
using QuantizedArrays
v = collect(1:10);
qv = QuantizedArray(v, k=3)  # 3 quantization levels
# 10-element QuantizedArray{UInt8,Int64,1}:
#  2
#  2
#  2
#  2
#  6
#  6
#  6
#  9
#  9
#  9

m = reshape(collect(1:24), (6,4))
qm = QuantizedArray(m, k = 2, m = 3)  # 2 vectors, 3 codebooks
# 6×4 QuantizedArray{UInt8,Int64,2}:
#   1   1  13  13
#   2   2  14  14
#   9   9  15  15
#  10  10  16  16
#  11  11  17  17
#  12  12  18  18
```

## Features
To keep track with the latest features, please consult [NEWS.md](https://github.com/zgornel/QuantizedArrays.jl/blob/master/NEWS.md) and the [documentation](https://zgornel.github.io/QuantizedArrays.jl/dev).


## License

The code has an MIT license and therefore it is free.


## Reporting Bugs

This is work in progress and bugs may still be present...¯\\_(ツ)_/¯ Do not worry, just [open an issue](https://github.com/zgornel/QuantizedArrays.jl/issues/new) to report a bug or request a feature.


## References

[1] The package [Rayuela.jl](https://github.com/una-dinosauria/Rayuela.jl) provides an extensive list of quantization methods implementations as well as good references pertinent to the state-of-the-art.

[2] [Quantization in signal processing (Wikipedia)](https://en.wikipedia.org/wiki/Quantization_(signal_processing))

[3] [Vector quantization (Wikipedia)](https://en.wikipedia.org/wiki/Vector_quantization)
