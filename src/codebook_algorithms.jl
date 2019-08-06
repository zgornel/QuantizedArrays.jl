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
"""
function opq_codebooks(X::AbstractMatrix{T}, k::Int, m::Int;
                       distance::Distances.PreMetric=DEFAULT_DISTANCE,
                       maxiter::Int=DEFAULT_OPQ_MAXITER) where {T}
    # Initialize R, rotated data
    nrows, ncols = size(X)
    X̂ = zeros(T, size(X))
    R = diagm(0 => ones(T, nrows))
    Xr = R * X

    # Initialize codebooks, codes
    cweights   = zeros(Int, k)
    costs      = zeros(ncols)
    counts     = zeros(Int, k)
    unused     = Int[]
    to_update_z = zeros(Bool, k)
    codes = fill(zeros(Int, ncols), m)
    cbooks = sampling_codebooks(X, k, m)
    @inbounds for i in 1:m
        rr = rowrange(nrows, m, i)
        dists = Distances.pairwise(distance, cbooks[i], Xr[rr, :], dims=2)
        Clustering.update_assignments!(dists, true, codes[i], costs, counts, to_update_z, unused)
        X̂[rr,:] .= cbooks[i][:, codes[i]]
    end

    # Run optimization
    to_update = fill(ones(Bool, k), m)
    for it = 1:maxiter
        # Update R using orthogonal Procustes closed form solution
        # and update rotated data matrix X̂
        U, _, V = svd(X * X̂', full=false)
        R = V * U'
        Xr = R * X
        @inbounds for i in 1:m
            rr = rowrange(nrows, m, i)
            # Update subspace cluster centers
            Clustering.update_centers!(Xr[rr, :], nothing, codes[i], to_update[i], cbooks[i], cweights)
            # Update subspace data assignments
            dists = Distances.pairwise(distance, cbooks[i], Xr[rr, :], dims=2)
            Clustering.update_assignments!(dists, false, codes[i], costs, counts, to_update[i], unused)
            X̂[rr, :] .= cbooks[i][:, codes[i]]
        end
        @debug "OPQ: Iteration $it, error = $(sum((R'*X̂ .- X).^2)./ncols)"
    end
    return cbooks, R
end


"""
    rvq_codebooks(X, k, m [;distance=DEFAULT_DISTANCE, maxiter=DEFAULT_RVQ_MAXITER])

Build `m` codebooks for input matrix `X` using `k` layers of residual vector cluster centers,
obtained using k-means clustering, the distance `distance` and `maxiter` iterations.

# References:
  * [Chen et al. 2010](https://www.mdpi.com/1424-8220/10/12/11259/htm)
"""
function rvq_codebooks(X::AbstractMatrix{T}, k::Int, m::Int;
                       distance::Distances.PreMetric=DEFAULT_DISTANCE,
                       maxiter::Int=DEFAULT_RVQ_MAXITER) where {T}
    nrows = size(X, 1)
    cbooks = Vector{Matrix{T}}(undef, m)
    Xr = copy(X)
    @inbounds for i in 1:m
        model = kmeans(Xr, k,
                       maxiter=maxiter,
                       distance=distance,
                       init=:kmpp,
                       display=:none)
        cbooks[i] = model.centers
        codes = model.assignments
        Xr .-= cbooks[i][:, codes]
    end
    return cbooks
end
