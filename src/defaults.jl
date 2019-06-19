# Default parameter values
const DEFAULT_K = 256
const DEFAULT_M = 1
const DEFAULT_METHOD = :sample
const DEFAULT_DISTANCE = Distances.SqEuclidean()
const MAX_BITS = 128
const BITS_TO_TYPE = Dict(8 => UInt8,
                          16 => UInt16,
                          32 => UInt32,
                          64 => UInt64,
                          128 => UInt128)
const TYPE_TO_BITS = Dict(t=>b for (b,t) in BITS_TO_TYPE)
const DEFAULT_PQ_MAXITER = 25
const DEFAULT_OPQ_MAXITER = 25
