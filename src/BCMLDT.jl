"""
    BCMLDT

A link-based solution method for the stochastic user equilibrium (SUE) under the
bounded choice model with local detour threshold (BCM-LDT), originally proposed by
Rasmussen et al. (2024).

Reference: Rasmussen, T.K., Duncan, L.C., Watling, D.P., Nielsen, O.A. (2024).
Local detouredness: A new phenomenon for modelling route choice and traffic
assignment. Transportation Research Part B 190, 103052.
"""
module BCMLDT

using LinearAlgebra: norm, tr, I
using Printf: @printf
using Downloads

export Network, ODPairs, BPRParams, BCMLDTParams
export num_links, num_ods, link_cost, link_cost!
export min_cost_table, min_cost_table_mt
export DFSWorkspace, pruned_dfs!
export bcmldt_loading!, bcmldt_loading, bcmldt_loading_mt!
export solve_bcmldt_sue, solve_bcmldt_sue_mt
export read_tntp_net, read_tntp_trips, read_tntp_ue, download_tntp

# engine (each file ends with a multi-threaded version, with the suffix _mt)
include("types.jl")              # inputs: network, OD pairs, parameters
include("shortest.jl")           # minimum cost matrix D
include("dfs.jl")                # pruned DFS (Algorithm 1.1)
include("loading.jl")            # BCM-LDT loading (Algorithm 1)
include("solver.jl")             # Anderson acceleration (Algorithm 2)

# utilities
include("util_tntp.jl")          # reading and downloading TNTP files

end
