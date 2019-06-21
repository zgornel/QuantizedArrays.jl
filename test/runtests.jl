module TestQuantizedArrays

using Test
using InteractiveUtils
using Distances
using QuantizedArrays

include("codebook.jl")
include("quantizer.jl")
include("interface.jl")
include("utils.jl")

end
