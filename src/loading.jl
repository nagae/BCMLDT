# =====================================================================
#  BCM-LDT loading L(t) (Algorithm 1)
# =====================================================================

"""
    bcmldt_loading!(x, net, od, t, D, prm) -> num_paths

BCM-LDT loading ``x = L(t)`` under link costs `t` and minimum cost matrix
`D = min_cost_table(net, t)` (Algorithm 1). The link flows are written into `x`,
and the number of paths with positive weight for each OD pair (the choice set size)
is returned.
For each OD pair ``w``, the pruned DFS gives ``Π^w_a`` and ``Π^w`` without storing
paths, and ``x_a`` is increased by ``q_w Π^w_a / Π^w`` (Eq. (7)).
The minimum cost of every OD pair must be positive and finite; when `γ` is finite,
all link costs must be positive so that the local detouredness is defined.
"""
function bcmldt_loading!(x::AbstractVector{Float64}, net::Network, od::ODPairs,
                         t::AbstractVector, D::Matrix{Float64}, prm::BCMLDTParams)
    @assert all(w -> 0 < D[od.orig[w], od.dest[w]] < Inf, 1:num_ods(od)) "minimum costs of OD pairs must be positive and finite"

    ws = DFSWorkspace(net.num_nodes)
    Πa = zeros(num_links(net))
    num_paths = zeros(Int, num_ods(od))
    fill!(x, 0.0)                                          # Step 1
    for w in 1:num_ods(od)                                 # Step 2
        fill!(Πa, 0.0)
        Π, num_paths[w] = pruned_dfs!(Πa, ws, net, t, D, prm, od.orig[w], od.dest[w])
        c = od.q[w] / Π
        @inbounds for a in eachindex(x)
            x[a] += c * Πa[a]                              # Eq. (7)
        end
    end
    return num_paths                                       # Step 3
end

"""
    bcmldt_loading(net, od, t, D, prm) -> (x, num_paths)

Same as `bcmldt_loading!`, but returns new link flows `x`.
"""
function bcmldt_loading(net::Network, od::ODPairs, t::AbstractVector,
                        D::Matrix{Float64}, prm::BCMLDTParams)
    x = zeros(num_links(net))
    num_paths = bcmldt_loading!(x, net, od, t, D, prm)
    return x, num_paths
end

# =====================================================================
#  Multi-threaded version of bcmldt_loading!
#
#  The OD pairs are processed in batches. Within a batch, the pruned DFS for each
#  OD pair is run in parallel, and its contribution q_w Π^w_a / Π^w is stored in
#  a column of a matrix. The columns are then added to x in the order of the OD
#  pairs, so that the result does not depend on the number of threads and is
#  identical to that of bcmldt_loading!.
# =====================================================================

"""
    bcmldt_loading_mt!(x, net, od, t, D, prm; batch_size = 1000) -> num_paths

Multi-threaded version of `bcmldt_loading!`, with the same result. The OD pairs are
processed in batches of `batch_size`, which bounds the memory for the contributions
of the OD pairs to `num_links(net) × batch_size` numbers.
Start Julia with multiple threads (e.g. `julia -t 8`) to make use of it.
"""
function bcmldt_loading_mt!(x::AbstractVector{Float64}, net::Network, od::ODPairs,
                            t::AbstractVector, D::Matrix{Float64}, prm::BCMLDTParams;
                            batch_size::Int = 1000)
    @assert all(w -> 0 < D[od.orig[w], od.dest[w]] < Inf, 1:num_ods(od)) "minimum costs of OD pairs must be positive and finite"

    contrib = zeros(num_links(net), min(batch_size, num_ods(od)))   # contributions of OD pairs
    num_paths = zeros(Int, num_ods(od))
    fill!(x, 0.0)
    for batch in Iterators.partition(1:num_ods(od), batch_size)
        # pruned DFS for the OD pairs in the batch, in parallel
        next = Threads.Atomic{Int}(1)
        @sync for _ in 1:Threads.nthreads()
            # the workspace is allocated within each task, so that the workspaces of
            # different tasks do not share cache lines (false sharing)
            Threads.@spawn begin
                ws = DFSWorkspace(net.num_nodes)
                while true
                    k = Threads.atomic_add!(next, 1)        # k-th OD pair in the batch
                    k > length(batch) && break
                    w = batch[k]
                    Πa = view(contrib, :, k)
                    fill!(Πa, 0.0)
                    Π, num_paths[w] = pruned_dfs!(Πa, ws, net, t, D, prm, od.orig[w], od.dest[w])
                    c = od.q[w] / Π
                    @inbounds for a in eachindex(Πa)
                        Πa[a] = c * Πa[a]
                    end
                end
            end
        end
        # add the contributions to x in the order of the OD pairs
        for k in eachindex(batch)
            @inbounds for a in eachindex(x)
                x[a] += contrib[a, k]
            end
        end
    end
    return num_paths
end
