# Codebook building algorithms
# ----------------------------
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
