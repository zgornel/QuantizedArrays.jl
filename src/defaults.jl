# Default parameter values
const DEFAULT_K = 256
const DEFAULT_M = 1
const DEFAULT_METHOD = :sample
const DEFAULT_DISTANCE = Distances.SqEuclidean()
const MAX_BITS = 128
const BITS_TO_TYPE = Dict(round(Int, log2(typemax(typ))) => typ
                          for typ in subtypes(Unsigned))
const TYPE_TO_BITS = Dict(t=>b for (b,t) in BITS_TO_TYPE)
const DEFAULT_PQ_MAXITER = 25
const DEFAULT_OPQ_MAXITER = 25
const DEFAULT_RVQ_MAXITER = 25
