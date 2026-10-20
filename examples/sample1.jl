# =====================================================================
#  Sample 1: choice sets of OD pair (4, 7) of the Anaheim network
#
#  Reproduces the comparison with Rasmussen et al. (2024) in Section 4: with link
#  travel times fixed at the UE values provided in the data, the numbers of paths
#  with positive weight for OD pair (4, 7) are 1, 1, 6, 33, 600, and 3,725 for
#  γ = 0.1, 0.2, 0.4, 0.6, 0.8, and 1.0.
#
#  Run from the root of the package:  julia --project=. examples/sample1.jl
# =====================================================================
using BCMLDT

dir = download_tntp("Anaheim", joinpath(@__DIR__, "..", "data"))
(; net, len) = read_tntp_net(joinpath(dir, "Anaheim_net.tntp"))
ue = read_tntp_ue(joinpath(dir, "Anaheim_flow.tntp"), net)

# link costs: UE travel time (min) + 0.5 × length (km); lengths are in feet and
# converted with 0.3/1000 km/ft, as in the implementation of Rasmussen et al. (2024)
t = ue.cost .+ 0.5 .* (0.3e-3 .* len)
D = min_cost_table(net, t)

r, s = 4, 7
ws = DFSWorkspace(net.num_nodes)
println("Number of paths with positive weight for OD pair ($(r), $(s)) (τ = 1.6)")
for γ in (0.1, 0.2, 0.4, 0.6, 0.8, 1.0)
    prm = BCMLDTParams(θ1 = 0.2, θ2 = 0.2, τ = 1.6, γ = γ)
    _, num_paths = pruned_dfs!(zeros(num_links(net)), ws, net, t, D, prm, r, s)
    println("  γ = $(γ): $(num_paths)")
end
