```@meta
CurrentModule=QuantizedArrays
```

# Introduction

`QuantizedArrays` is a package for array quantization and compression. The basic principle is that arrays can be represented through the concatenation of shorter vectors obtained by either sampling or clustering of subspaces in the original data. This effectively compresses the original array with a loss in quality.

## Implemented features
 - basic `AbstractArray` interface for vectors, matrices
 - quantization types implemented:
     - random sampling
     - k-mean clustering ([PQ](https://lear.inrialpes.fr/pubs/2011/JDS11/jegou_searching_with_quantization.pdf))
     - 'cartesian' k-means clustering ([OPQ](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/opq_tr.pdf))
     - residual quantization ([RVQ](https://www.mdpi.com/1424-8220/10/12/11259/htm))

## Installation

Installation can be performed from either inside or outside Julia.

### Git cloning
```
$ git clone https://github.com/zgornel/QuantizedArrays.jl
```

### Julia REPL
The package can be installed from inside Julia with:
```
using Pkg
Pkg.add("QuantizedArrays")
```
or
```
Pkg.add(PackageSpec(url="https://github.com/zgornel/QuantizedArrays.jl", rev="master"))
```
for the latest `master` branch.
