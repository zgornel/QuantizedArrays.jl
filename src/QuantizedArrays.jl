__precompile__(true)

module QuantizedArrays

export QuantizedVector,
       QuantizedMatrix,
       QuantizedArray,
       build_quantizer,
       quantize,
       quantizer,
       codebook

include("defaults.jl")
include("quantize.jl")
include("interface.jl")

end # module
