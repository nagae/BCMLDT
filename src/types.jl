# =====================================================================
#  Types: network, OD pairs, parameters of link cost functions and of BCM-LDT
# =====================================================================

"""
    Network(num_nodes, tail, head; first_thru_node = 1)

Directed graph ``G = (N, A)``. Nodes are numbered `1:num_nodes`, and link `a` goes from
`tail[a]` to `head[a]`; links are numbered in the order given. Nodes with numbers
smaller than `first_thru_node` (the OD nodes, or zones) are used only as origins and
destinations and are never passed through.

The links leaving node `i` are `out_link[out_ptr[i]:out_ptr[i+1]-1]`, in the order given.
"""
struct Network
    num_nodes::Int
    first_thru_node::Int
    tail::Vector{Int}
    head::Vector{Int}
    # links grouped by their tail node; those leaving node i are out_link[out_ptr[i]:out_ptr[i+1]-1]
    out_link::Vector{Int}
    out_ptr::Vector{Int}
    # head of each out_link: out_head[e] == head[out_link[e]]
    # (a copy kept contiguous for speed in the DFS and Dijkstra's algorithm)
    out_head::Vector{Int}
end

function Network(num_nodes::Integer, tail::AbstractVector{<:Integer},
                 head::AbstractVector{<:Integer}; first_thru_node::Integer = 1)
    @assert all(i -> 1 <= i <= num_nodes, head) "node numbers must be in 1:num_nodes"

    out_link = sortperm(tail; alg = MergeSort)     # stable: keeps the link order within a node
    out_ptr = zeros(Int, num_nodes + 1)
    for i in tail
        out_ptr[i+1] += 1
    end
    out_ptr[1] = 1
    cumsum!(out_ptr, out_ptr)
    out_head = head[out_link]

    # default constructor (all fields)
    return Network(num_nodes, first_thru_node, collect(Int, tail), collect(Int, head),
                   out_link, out_ptr, out_head)
end

"Number of links ``|A|``."
num_links(net::Network) = length(net.tail)

"Whether node `i` is an OD node (zone), which is never passed through."
is_zone(net::Network, i::Integer) = i < net.first_thru_node

"""
    ODPairs(orig, dest, q)
    ODPairs(Q::AbstractMatrix)

OD pairs ``w = (r, s)`` with demand ``q_w``. From a matrix `Q[r, s]`, OD pairs with
zero demand and with `r == s` are dropped, and the pairs are ordered by origin.
"""
struct ODPairs
    orig::Vector{Int}
    dest::Vector{Int}
    q::Vector{Float64}
end

function ODPairs(Q::AbstractMatrix)
    orig, dest, q = Int[], Int[], Float64[]
    for r in axes(Q, 1), s in axes(Q, 2)          # s is the inner loop: OD pairs are grouped by origin
        (r == s || Q[r, s] == 0) && continue
        push!(orig, r); push!(dest, s); push!(q, Q[r, s])
    end
    return ODPairs(orig, dest, q)
end

"Number of OD pairs ``|W|``."
num_ods(od::ODPairs) = length(od.q)

"""
    BPRParams(t0, μ, α, β, t_add)

Parameters of the link cost function of the Bureau of Public Roads (BPR) type with a
flow-independent additional cost:
``t_a(x_a) = t0_a \\{1 + α_a ([x_a]_+ / μ_a)^{β_a}\\} + t_add_a``,
where `t0` is the free flow time, `μ` the capacity, `α` and `β` the BPR coefficient
and exponent (often 0.15 and 4), and `t_add` an additional cost (e.g. a
distance-based cost). Link flows are replaced by their positive parts ``[x_a]_+``,
so that the cost is defined for any ``x`` (see the footnote to Step 1 of Algorithm 2).
"""
struct BPRParams
    t0::Vector{Float64}
    μ::Vector{Float64}
    α::Vector{Float64}
    β::Vector{Float64}
    t_add::Vector{Float64}
end

"Compute link costs ``t(x)`` into `t`."
function link_cost!(t::AbstractVector, cost_params::BPRParams, x::AbstractVector)
    (; t0, μ, α, β, t_add) = cost_params
    @inbounds for a in eachindex(t)
        t[a] = t0[a] * (1 + α[a] * (max(x[a], 0.0) / μ[a])^β[a]) + t_add[a]
    end
    return t
end

"Link costs ``t(x)``."
link_cost(cost_params::BPRParams, x::AbstractVector) = link_cost!(similar(cost_params.t0), cost_params, x)

"""
    BCMLDTParams(; θ1, θ2, τ, γ)

Parameters of BCM-LDT: scaling parameters for travel cost `θ1` and local
detouredness `θ2`, global cost bound `τ > 1`, and local detour threshold `γ > 0`.
`θ2 = 0` with `γ = Inf` gives the cost-based BCM.
"""
struct BCMLDTParams
    θ1::Float64
    θ2::Float64
    τ::Float64
    γ::Float64
    function BCMLDTParams(θ1, θ2, τ, γ)
        θ1 > 0 || throw(ArgumentError("θ1 must be positive"))
        θ2 >= 0 || throw(ArgumentError("θ2 must be nonnegative"))
        τ > 1 || throw(ArgumentError("τ must be greater than 1"))
        γ > 0 || throw(ArgumentError("γ must be positive"))
        isfinite(γ) && θ2 == 0 && throw(ArgumentError("θ2 must be positive when γ is finite"))
        return new(θ1, θ2, τ, γ)
    end
end
BCMLDTParams(; θ1, θ2, τ, γ) = BCMLDTParams(θ1, θ2, τ, γ)
