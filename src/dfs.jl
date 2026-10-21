# =====================================================================
#  Pruned depth-first search for one OD pair (Algorithm 1.1)
# =====================================================================

# Working arrays for the DFS of one OD pair.
mutable struct DFSWorkspace
    node::Vector{Int}          # prefix (i_1, ..., i_n)
    cum::Vector{Float64}       # cum[l] = C(r → i_l)
    on_path::Vector{Bool}      # on_path[i]: whether node i is on the prefix (C1)
    Π::Float64                 # total weight of the paths found, summed in the order found
                               # (equal to the value returned for the prefix (r) up to rounding)
    num_paths::Int             # number of paths found with positive weight
end
DFSWorkspace(num_nodes::Int) =
    DFSWorkspace(Vector{Int}(undef, num_nodes), Vector{Float64}(undef, num_nodes),
                 zeros(Bool, num_nodes), 0.0, 0)

# exp(x) overflows Float64 for x > log(floatmax(Float64)) ≈ 709.78; keep a margin
const _EXP_ARG_MAX = 700.0

# Weight π_k of a path with cost C and local detouredness φ (Eq. (5)), divided by
# a constant common to all paths of the OD pair to avoid overflow and underflow
# (this does not change the choice probabilities P_k = π_k / Π^w):
#   (e^{-θ1 (C - τ d_w)} - 1) / e^{θ1 (τ - 1) d_w}  = e^{θ1 (d_w - C)} - e^{θ1 (d_w - τ d_w)},
#   (e^{-θ2 (φ - γ)} - 1)                           (divided by e^{θ2 γ} if θ2 γ ≥ _EXP_ARG_MAX).
# With γ = Inf (cost-based BCM), the second factor is omitted.
@inline function _path_weight(C, φ, d_w, g_bound, prm::BCMLDTParams)
    (; θ1, θ2, γ) = prm
    u = exp(θ1 * (d_w - C)) - exp(θ1 * (d_w - g_bound))
    isfinite(γ) || return u
    # e^{θ2 (γ - φ)} - 1 by expm1, which avoids cancellation when φ ≈ γ; if θ2 γ is large
    # enough for exp to overflow, the factor divided by e^{θ2 γ} instead
    v = θ2 * γ < _EXP_ARG_MAX ? expm1(θ2 * (γ - φ)) : exp(-θ2 * φ) - exp(-θ2 * γ)
    return u * v
end

# PrunedDFS of Algorithm 1.1. The prefix (i_1, ..., i_n) is ws.node[1:n], with cost
# C = C(r → i_n) and local detouredness φ = φ(r → i_n). Adds to Πa[a] the total weight
# of the paths using link a, and returns the total weight of the paths sharing the prefix.
function _pruned_dfs!(Πa::AbstractVector{Float64}, ws::DFSWorkspace, net::Network,
                      t::AbstractVector, D::Matrix{Float64}, prm::BCMLDTParams,
                      s::Int, d_w::Float64, g_bound::Float64, n::Int, C::Float64, φ::Float64)
    i = ws.node[n]
    # Step 1: a path has been found
    if i == s
        w = _path_weight(C, φ, d_w, g_bound, prm)
        ws.Π += w
        w > 0 && (ws.num_paths += 1)                       # a path on the boundary has π_k = 0
        return w
    end
    # Step 2
    σ = 0.0
    # Step 3: extensions of the prefix by the links leaving i
    @inbounds for e in net.out_ptr[i]:net.out_ptr[i+1]-1
        j = net.out_head[e]; a = net.out_link[e]
        ws.on_path[j] && continue                          # (C1)
        is_zone(net, j) && j != s && continue              # OD nodes are not passed through
        C_next = C + t[a]                                      # Step 3.1: C_next = C(r → j)
        C_next + D[j, s] > g_bound && continue                 # (C2): C(r → j) + d(j, s) > τ d_w
        φ_next = φ                                             # Step 3.3 and (C3)
        if isfinite(prm.γ)
            pruned = false
            for l in n:-1:1                                # nearest to j first: prunes earlier
                d = D[ws.node[l], j]
                ρ = (C_next - ws.cum[l] - d) / d               # (C(i_l → j) - d(i_l, j)) / d(i_l, j)
                if ρ > prm.γ
                    pruned = true; break
                end
                ρ > φ_next && (φ_next = ρ)
            end
            pruned && continue
        end
        ws.node[n+1] = j; ws.cum[n+1] = C_next; ws.on_path[j] = true   # Step 3.4
        σ_child = _pruned_dfs!(Πa, ws, net, t, D, prm, s, d_w, g_bound, n + 1, C_next, φ_next)
        ws.on_path[j] = false
        Πa[a] += σ_child                                        # Step 3.5
        σ += σ_child
    end
    # Step 4
    return σ
end

"""
    pruned_dfs!(Πa, ws, net, t, D, prm, r, s) -> (Π, num_paths)

Pruned depth-first search of Rasmussen et al. (2024) for OD pair ``w = (r, s)``
under link costs `t` and minimum cost matrix `D` (Algorithm 1.1, called with the
prefix ``(r)``). Adds to `Πa[a]` the total weight ``Π^w_a`` of the paths using link
`a`, and returns the total weight ``Π^w`` of all paths and the number of paths with
positive weight (the choice set size). Paths on the boundary of the bounds
(``C_k = τ d_w`` or ``φ_k = γ``) may also be found, but they have zero weight.
Paths are enumerated but not stored. The weights are divided by a constant common to
all paths of the OD pair, which does not change the choice probabilities.
"""
function pruned_dfs!(Πa::AbstractVector{Float64}, ws::DFSWorkspace, net::Network,
                     t::AbstractVector, D::Matrix{Float64}, prm::BCMLDTParams,
                     r::Int, s::Int)
    d_w = D[r, s]
    ws.node[1] = r; ws.cum[1] = 0.0; ws.on_path[r] = true
    ws.Π = 0.0; ws.num_paths = 0
    _pruned_dfs!(Πa, ws, net, t, D, prm, s, d_w, prm.τ * d_w, 1, 0.0, 0.0)   # g_bound = τ d_w
    ws.on_path[r] = false
    return ws.Π, ws.num_paths
end
