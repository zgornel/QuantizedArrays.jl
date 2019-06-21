"""
    sampling_codebooks(X, k, m)

Build `m` codebooks for input matrix `X` by sampling `k` prototype vectors.
"""
function sampling_codebooks(X::AbstractMatrix{T}, k::Int, m::Int) where {T}
    nrows, ncols = size(X)
    cbooks = Vector{Matrix{T}}(undef, m)
    @inbounds for i in 1:m
        rr = rowrange(nrows, m, i)
        vectors = Matrix{T}(X[rr, sample(1:ncols, k, replace=false)])
        cbooks[i] = vectors
    end
    return cbooks
end


"""
    pq_codebooks(X, k, m [;distance=DEFAULT_DISTANCE, maxiter=DEFAULT_PQ_MAXITER])

Build `m` codebooks for input matrix `X` using `k` sub-space centers obtained
using k-means clustering, the distance `distance` and `maxiter` iterations.

# References:
  * [Jègou et al. 2011](https://lear.inrialpes.fr/pubs/2011/JDS11/jegou_searching_with_quantization.pdf)
#
"""
function pq_codebooks(X::AbstractMatrix{T}, k::Int, m::Int;
                      distance::Distances.PreMetric=DEFAULT_DISTANCE,
                      maxiter::Int=DEFAULT_PQ_MAXITER) where {T}
    nrows = size(X, 1)
    cbooks = Vector{Matrix{T}}(undef, m)
    @inbounds for i in 1:m
        rr = rowrange(nrows, m, i)
        model = kmeans(X[rr, :], k,
                       maxiter=maxiter,
                       distance=distance,
                       init=:kmpp,
                       display=:none)
        cbooks[i] = model.centers
    end
    return cbooks
end


"""
    opq_codebooks(X, k, m [;distance=DEFAULT_DISTANCE, maxiter=DEFAULT_PQ_MAXITER])

Build `m` codebooks for input matrix `X` using `k` sub-space centers obtained
using 'cartesian' k-means clustering, the distance `distance` and `maxiter` iterations.

# References:
  * [Ge et al. 2014](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/opq_tr.pdf)
#
"""
function opq_codebooks(X::AbstractMatrix{T}, k::Int, m::Int;
                       distance::Distances.PreMetric=DEFAULT_DISTANCE,
                       maxiter::Int=DEFAULT_OPQ_MAXITER) where {T}
    # Initialize R, rotated data
    nrows, ncols = size(X)
    R = diagm(0 => ones(T, nrows))
    X̂ = R * X

    # Initialize codebooks, codes
    cweights   = zeros(Int, k)
    costs      = zeros(ncols)
    counts     = zeros(Int, k)
    unused     = Int[]
    to_update = zeros(Bool, k)
    codes = fill(zeros(Int, ncols), m)
    cbooks = pq_codebooks(X, k, m, distance=distance, maxiter=10)
    @inbounds for i in 1:m
        rr = rowrange(nrows, m, i)
        dists = Distances.pairwise(distance, cbooks[i], X̂[rr, :], dims=2)
        Clustering.update_assignments!(dists, false, codes[i], costs, counts, to_update, unused)
        X̂[rr,:] .= cbooks[i][:, codes[i]]
    end

    # Run optimization
    to_update = ones(Bool, k)
    for _ = 1:maxiter
        # Update R using orthogonal Procustes closed form solution
        # and update rotated data matrix X̂
        U, _, V = svd(X * X̂', full=false)
        R = V * U'
        X̂ = R * X
        @inbounds for i in 1:m
            rr = rowrange(nrows, m, i)
            # Update subspace cluster centers
            Clustering.update_centers!(X̂[rr, :], nothing, codes[i], to_update, cbooks[i], cweights)
            # Update subspace data assignments
            dists = Distances.pairwise(distance, cbooks[i], X̂[rr, :], dims=2)
            Clustering.update_assignments!(dists, false, codes[i], costs, counts, to_update, unused)
            X̂[rr, :] .= cbooks[i][:, codes[i]]
        end
    end
    return cbooks
end
