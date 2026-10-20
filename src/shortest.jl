# =====================================================================
#  Minimum cost matrix D = [d(i, j)] (Algorithm 2, Step 1)
#
#  Dijkstra's algorithm is implemented here instead of using Graphs.jl, because
#  OD nodes (zones) must not be passed through, which `dijkstra_shortest_paths`
#  of Graphs.jl does not support.
# =====================================================================

# Dijkstra's algorithm from `src`. The priority queue is a binary heap (hk: keys =
# tentative costs, hv: nodes). When the tentative cost of a node decreases, a new
# entry is pushed instead of updating the old one; outdated entries are skipped when
# popped ("lazy deletion").
# OD nodes (zones) are not expanded, i.e., they are never passed through;
# this includes `src` itself when it is a zone (see `min_cost_table`).
function _dijkstra!(dist::Vector{Float64}, hk::Vector{Float64}, hv::Vector{Int},
                    net::Network, t::AbstractVector, src::Int)
    fill!(dist, Inf)
    dist[src] = 0.0
    n = 1; hk[1] = 0.0; hv[1] = src
    @inbounds while n > 0
        du, u = hk[1], hv[1]                       # pop the minimum
        hk[1], hv[1] = hk[n], hv[n]; n -= 1
        i = 1
        while true                                 # sift down
            l = 2i; r = l + 1; m = i
            l <= n && hk[l] < hk[m] && (m = l)
            r <= n && hk[r] < hk[m] && (m = r)
            m == i && break
            hk[i], hk[m] = hk[m], hk[i]; hv[i], hv[m] = hv[m], hv[i]
            i = m
        end
        du > dist[u] && continue                   # outdated entry
        is_zone(net, u) && continue                # zones are not passed through
        for e in net.out_ptr[u]:net.out_ptr[u+1]-1
            v = net.out_head[e]
            nd = du + t[net.out_link[e]]
            if nd < dist[v]
                dist[v] = nd
                n += 1                             # push
                if n > length(hk)
                    resize!(hk, 2n); resize!(hv, 2n)
                end
                hk[n] = nd; hv[n] = v
                j = n
                while j > 1                        # sift up
                    pj = j >> 1
                    hk[pj] <= hk[j] && break
                    hk[pj], hk[j] = hk[j], hk[pj]; hv[pj], hv[j] = hv[j], hv[pj]
                    j = pj
                end
            end
        end
    end
    return dist
end

"""
    min_cost_table(net, t) -> D

Minimum cost matrix: `D[i, j]` is the minimum cost ``d(i, j)`` from node `i` to
node `j` under link costs `t`, where OD nodes (zones) are not passed through.
Rows for non-zone nodes are computed by Dijkstra's algorithm. A zone `r` is left
only through its outgoing links, so its row is
``d(r, k) = \\min_{a = (r, j)} \\{t_a + d(j, k)\\}``.
Link costs must be nonnegative. See `min_cost_table_mt` for a multi-threaded version.
"""
function min_cost_table(net::Network, t::AbstractVector)
    nn = net.num_nodes
    D = Matrix{Float64}(undef, nn, nn)
    dist = Vector{Float64}(undef, nn)
    # heap for Dijkstra's algorithm: at most one entry per link plus the source
    hk = Vector{Float64}(undef, num_links(net) + 1)
    hv = Vector{Int}(undef, num_links(net) + 1)
    for i in 1:nn
        _dijkstra!(dist, hk, hv, net, t, i)
        D[i, :] .= dist
    end
    _zone_rows!(D, net, t)
    return D
end

# Rows of D for zones: d(r, k) = min_{a = (r, j)} {t_a + d(j, k)}.
function _zone_rows!(D::Matrix{Float64}, net::Network, t::AbstractVector)
    nn = net.num_nodes
    for r in 1:net.first_thru_node-1
        best = fill(Inf, nn); best[r] = 0.0
        for e in net.out_ptr[r]:net.out_ptr[r+1]-1
            j = net.out_head[e]; ta = t[net.out_link[e]]
            @inbounds for k in 1:nn
                v = ta + D[j, k]
                v < best[k] && (best[k] = v)
            end
        end
        D[r, :] .= best
    end
    return D
end

# =====================================================================
#  Multi-threaded version of min_cost_table
#
#  The rows of D are computed in parallel over nodes. Each row is computed
#  independently, so the result does not depend on the number of threads and is
#  identical to that of min_cost_table.
# =====================================================================

"""
    min_cost_table_mt(net, t) -> D

Multi-threaded version of `min_cost_table`, with the same result.
Start Julia with multiple threads (e.g. `julia -t 8`) to make use of it.
"""
function min_cost_table_mt(net::Network, t::AbstractVector)
    nn = net.num_nodes
    D = Matrix{Float64}(undef, nn, nn)
    Threads.@threads for i in 1:nn
        dist = Vector{Float64}(undef, nn)
        # heap for Dijkstra's algorithm: at most one entry per link plus the source
        hk = Vector{Float64}(undef, num_links(net) + 1)
        hv = Vector{Int}(undef, num_links(net) + 1)
        _dijkstra!(dist, hk, hv, net, t, i)
        @inbounds for j in 1:nn
            D[i, j] = dist[j]
        end
    end
    _zone_rows!(D, net, t)
    return D
end
