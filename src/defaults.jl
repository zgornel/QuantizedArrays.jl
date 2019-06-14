# Default parameter values
const DEFAULT_K = 256
const DEFAULT_M = 1
const DEFAULT_METHOD = :sampling
const DEFAULT_DISTANCE = Distances.Euclidean()
const MAX_BITS = 128
const BITS_TO_TYPE = Dict(8 => UInt8,
                          16 => UInt16,
                          32 => UInt32,
                          64 => UInt64,
                          128 => UInt128)