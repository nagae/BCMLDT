# =====================================================================
#  Sample 2: BCM-LDT SUE on the Anaheim network
#
#  Reproduces Table 1 (and the data of Fig. 1) in Section 4: BCM-LDT SUE is solved
#  by Algorithm 2 for γ = 0.6, 0.8, and 1.0. Also shows the number of paths of
#  OD pair (4, 7) at the equilibrium with γ = 0.8 (351, as reported by
#  Rasmussen et al., 2024).
#
#  Run from the root of the package:  julia --project=. examples/sample2.jl
#  With multiple threads (e.g. julia -t 8 --project=. examples/sample2.jl),
#  the multi-threaded version solve_bcmldt_sue_mt is used.
#  The computation takes a few minutes with a single thread.
# =====================================================================
using BCMLDT, Printf

dir = download_tntp("Anaheim", joinpath(@__DIR__, "..", "data"))
(; net, cost_params, len) = read_tntp_net(joinpath(dir, "Anaheim_net.tntp"))
od = read_tntp_trips(joinpath(dir, "Anaheim_trips.tntp"))

# link costs: BPR travel time (min) + 0.5 × length (km); lengths are in feet and
# converted with 0.3/1000 km/ft, as in the implementation of Rasmussen et al. (2024)
(; t0, μ, α, β) = cost_params
cost_params = BPRParams(t0, μ, α, β, 0.5 .* (0.3e-3 .* len))

solve = Threads.nthreads() == 1 ? solve_bcmldt_sue : solve_bcmldt_sue_mt
solve(net, od, cost_params, BCMLDTParams(θ1 = 0.2, θ2 = 0.2, τ = 1.6, γ = 0.6); max_iter = 1)  # compilation

println("threads = $(Threads.nthreads())")
println("γ     iterations   time (s)   share of Algorithm 1 (%)   choice set size")
results = Dict{Float64,Any}()
for γ in (0.6, 0.8, 1.0)
    prm = BCMLDTParams(θ1 = 0.2, θ2 = 0.2, τ = 1.6, γ = γ)
    res = solve(net, od, cost_params, prm; m = 6, ω = 0.2, ε = 1e-9)
    results[γ] = res
    @printf("%.1f   %10d   %8.1f   %24.1f   %15d\n", γ, length(res.residuals) - 1, res.time_total,
            100 * res.time_loading / res.time_total, sum(res.num_paths))
end

w47 = findfirst(w -> (od.orig[w], od.dest[w]) == (4, 7), 1:length(od.q))
println("\nNumber of paths of OD pair (4, 7) at the equilibrium with γ = 0.8: $(results[0.8].num_paths[w47])")

println("\nRelative residuals ‖e_n‖/‖x_n‖ (Fig. 1)")
for γ in (0.6, 0.8, 1.0)
    println("γ = $(γ): ", join((@sprintf("%.2e", v) for v in results[γ].residuals), " "))
end
