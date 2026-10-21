# BCMLDT

A Julia implementation of the link-based solution method for the stochastic user
equilibrium (SUE) under the bounded choice model with local detour threshold (BCM-LDT)
of Rasmussen et al. (2024), described in

> T. Nagae, A link-based algorithm for the stochastic user equilibrium under the bounded
> choice model with local detour threshold, arXiv:XXXX.XXXXX.

The SUE is solved as a fixed point in aggregate link flows, x = L(t(x)), where the
loading L enumerates the paths in each choice set by a pruned depth-first search
(Algorithm 1) and the fixed-point problem is solved by Anderson acceleration (Algorithm 2).

## Requirements

Julia 1.10 or later. Only standard libraries are used.

## Examples

Run from the root of the package. The data of the Anaheim network are downloaded into
`data/` on the first run.

```
julia --project=. examples/sample1.jl
julia --project=. examples/sample2.jl
julia -t 8 --project=. examples/sample2.jl    # multi-threaded version
```

| File | Reproduces |
|---|---|
| `examples/sample1.jl` | Section 4: numbers of paths of OD pair (4, 7) with link travel times fixed at the UE values (1, 1, 6, 33, 600, 3,725 for γ = 0.1, …, 1.0) |
| `examples/sample2.jl` | Section 4: Table 1, the residual histories in Fig. 1, and the number of paths of OD pair (4, 7) at the equilibrium with γ = 0.8 (351) |

Computation times depend on the machine; iteration counts and numbers of paths do not.
The multi-threaded functions (suffix `_mt`) return the same results as the
single-threaded ones.

## Usage

The following solves the BCM-LDT SUE of the Anaheim network with γ = 0.6 in Table 1
(about 20 s with a single thread).

```julia
using BCMLDT

dir = download_tntp("Anaheim", "data")    # into ./data; files already there are not downloaded again
(; net, cost_params, len) = read_tntp_net(joinpath(dir, "Anaheim_net.tntp"))
od = read_tntp_trips(joinpath(dir, "Anaheim_trips.tntp"))

# link costs: BPR travel time (min) + 0.5 × length (km; lengths are in feet)
(; t0, μ, α, β) = cost_params
cost_params = BPRParams(t0, μ, α, β, 0.5 .* (0.3e-3 .* len))

prm = BCMLDTParams(θ1 = 0.2, θ2 = 0.2, τ = 1.6, γ = 0.6)
res = solve_bcmldt_sue(net, od, cost_params, prm)    # or solve_bcmldt_sue_mt
res.x                # equilibrium link flows
res.num_paths        # number of paths with positive weight for each OD pair
res.residuals        # relative residuals ‖e_n‖/‖x_n‖
```

The correspondence between the paper and the source files is

| Paper | Source |
|---|---|
| Algorithm 1.1 (pruned DFS) | `src/dfs.jl` |
| Algorithm 1 (loading), Eq. (7) | `src/loading.jl` |
| Algorithm 2 (Anderson acceleration) | `src/solver.jl` |
| minimum costs d(i, s) used in (C2) | `src/shortest.jl` |

See the docstrings (e.g. `?solve_bcmldt_sue`) for the arguments.

## Data

The networks are taken from the following repository, fixed at commit
`375e0da93858c547230c5cf9ea8a96de4ccff29e` (see `download_tntp`):

Transportation Networks for Research Core Team. Transportation Networks for Research.
https://github.com/bstabler/TransportationNetworks. Accessed October 20, 2026.

## Use of AI

The code was written with the assistance of Claude (Anthropic; Claude Opus 5.5 and
Claude Fable 5.1), and every function was reviewed by the author.

## License

MIT License (see `LICENSE`).

## Reference

Rasmussen, T.K., Duncan, L.C., Watling, D.P., Nielsen, O.A., 2024. Local detouredness:
A new phenomenon for modelling route choice and traffic assignment. Transportation
Research Part B 190, 103052.
