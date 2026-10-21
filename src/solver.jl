# =====================================================================
#  Solution of BCM-LDT SUE by fixed-point iteration with Anderson acceleration
#  (Algorithm 2)
# =====================================================================

# Parameters of the procedure used when the residual norm does not decrease in Step 4
# (see `_select_step_size`): the step size is halved at most _MAX_HALVING times, and if
# this fails, the stored differences are discarded and x moves by _ω_RESTART * e.
const _MAX_HALVING = 5
const _ω_RESTART = 0.02
# Tikhonov regularization of the least-squares problem for η (see Step 3)
const _TIKHONOV = 1e-10

"""
    solve_bcmldt_sue(net, od, cost_params, prm; m = 6, ω = 0.2, ε = 1e-9,
                     max_iter = 200, verbose = false)
        -> (x, num_paths, converged, residuals, time_total, time_loading)

Solve BCM-LDT SUE as the fixed-point problem in link flows
``x = F(x) = L(t(x))`` by Anderson acceleration (Algorithm 2), where ``t(x)`` is
given by `cost_params` and ``L`` is the BCM-LDT loading (`bcmldt_loading!`).
`m` is the number of stored differences, `ω` the damping parameter, and `ε` the
tolerance of the relative residual ``\\|F(x_n) - x_n\\| / \\|x_n\\|``.

Returns the equilibrium link flows `x` (``F(x_n)`` instead of ``x_n`` if ``x_n`` has
negative components), the number of paths for each OD pair at `x`
(the choice set sizes), whether the iteration converged, the relative residuals
for ``n = 0, 1, ...``, the total computation time, and the time spent in the
BCM-LDT loading (Algorithm 1).
See `solve_bcmldt_sue_mt` for a multi-threaded version.
"""
solve_bcmldt_sue(net::Network, od::ODPairs, cost_params::BPRParams, prm::BCMLDTParams; kwargs...) =
    _solve_bcmldt_sue(min_cost_table, bcmldt_loading!, net, od, cost_params, prm; kwargs...)

# Algorithm 2 with the functions for the minimum cost matrix and the loading given
# (single-threaded or multi-threaded versions).
function _solve_bcmldt_sue(min_cost_table_fn, loading_fn!, net::Network, od::ODPairs,
                           cost_params::BPRParams, prm::BCMLDTParams;
                           m::Int = 6, ω::Float64 = 0.2, ε::Float64 = 1e-9,
                           max_iter::Int = 200, verbose::Bool = false)
    t = zeros(num_links(net))
    time_loading = 0.0
    time_start = time()

    # Step 1 (Evaluation): F(x) = L(t(x)); link costs are evaluated at [x]_+.
    # Returns F(x) and the number of paths for each OD pair.
    function evaluate(x)
        link_cost!(t, cost_params, x)
        D = min_cost_table_fn(net, t)
        Fx = zeros(num_links(net))
        time_loading += @elapsed num_paths = loading_fn!(Fx, net, od, t, D, prm)
        return Fx, num_paths
    end

    # Step 0 (Initialization): x_0 = F(0), i.e., the loading under the free-flow
    # link costs t(0) (including the additional costs)
    x, _ = evaluate(zeros(num_links(net)))
    Fx, num_paths = evaluate(x)
    e = Fx .- x
    S = Vector{Vector{Float64}}()                          # s_i = x_{i+1} - x_i
    Y = Vector{Vector{Float64}}()                          # y_i = e_{i+1} - e_i
    residuals = Float64[]
    for n in 0:max_iter
        # Step 2 (Convergence test)
        push!(residuals, norm(e) / norm(x))
        verbose && @printf("n = %3d   relative residual = %.3e\n", n, residuals[end])
        if residuals[end] < ε
            any(<(0), x) && (x = x .+ e)                  # return F(x_n) = x_n + e_n if x_n has negative components
            return (x = x, num_paths = num_paths, converged = true, residuals = residuals,
                    time_total = time() - time_start, time_loading = time_loading)
        end
        n == max_iter && break

        # Step 3 (Direction finding): Eq. (8)
        #   d_n = ω e_n - Σ_i η_i (s_i + ω y_i),
        # where η = (η_i) solves the least-squares problem min_η ‖e_n - Σ_i η_i y_i‖.
        d = ω .* e
        if !isempty(S)
            Sm = reduce(hcat, S); Ym = reduce(hcat, Y)  # columns s_i and y_i (i = n-m, ..., n-1)
            # least squares by the normal equations (Y'Y + reg I) η = Y'e, with Tikhonov
            # regularization relative to the mean of ‖y_i‖² (= tr(Y'Y) / number of columns)
            YtY = Ym' * Ym
            reg = _TIKHONOV * tr(YtY) / size(YtY, 1)
            η = (YtY + reg * I) \ (Ym' * e)
            d .-= (Sm .+ ω .* Ym) * η
        end

        # Step 4 (Move): x_{n+1} = x_n + d_n
        x_next = x .+ d
        F_next, num_paths_next = evaluate(x_next)          # reused in Step 1 of the next iteration
        e_next = F_next .- x_next
        if norm(e_next) >= norm(e)                         # residual not reduced: adjust the step size
            x_next, F_next, e_next, num_paths_next, restarted = _select_step_size(evaluate, x, e, d)
            restarted && (empty!(S); empty!(Y))
        end

        # store the differences s_n and y_n, keeping the m most recent ones
        push!(S, x_next .- x); push!(Y, e_next .- e)
        if length(S) > m
            popfirst!(S); popfirst!(Y)
        end
        x, e, num_paths = x_next, e_next, num_paths_next
    end
    return (x = x, num_paths = num_paths, converged = false, residuals = residuals,
            time_total = time() - time_start, time_loading = time_loading)
end

# Step-size selection for Step 4 when the residual norm does not decrease with the step d:
# tries the step sizes 1/2, 1/4, ..., 1/2^_MAX_HALVING, and if none of them decreases
# the residual norm, moves by _ω_RESTART * e (restart). Returns the new point, F and
# the residual there, the number of paths, and whether it restarted.
function _select_step_size(evaluate, x, e, d)
    α = 1.0
    for _ in 1:_MAX_HALVING
        α /= 2
        x_next = x .+ α .* d
        F_next, num_paths_next = evaluate(x_next)
        e_next = F_next .- x_next
        norm(e_next) < norm(e) && return x_next, F_next, e_next, num_paths_next, false
    end
    x_next = x .+ _ω_RESTART .* e
    F_next, num_paths_next = evaluate(x_next)
    return x_next, F_next, F_next .- x_next, num_paths_next, true
end

# =====================================================================
#  Multi-threaded version of solve_bcmldt_sue
#
#  Only the minimum cost matrix and the BCM-LDT loading in Step 1, which take
#  almost all the computation time, are computed in parallel; the result is
#  identical to that of solve_bcmldt_sue.
# =====================================================================

"""
    solve_bcmldt_sue_mt(net, od, cost_params, prm; m = 6, ω = 0.2, ε = 1e-9,
                        max_iter = 200, verbose = false)

Multi-threaded version of `solve_bcmldt_sue`, with the same result. It uses
`min_cost_table_mt` and `bcmldt_loading_mt!` in Step 1.
Start Julia with multiple threads (e.g. `julia -t 8`) to make use of it.
"""
solve_bcmldt_sue_mt(net::Network, od::ODPairs, cost_params::BPRParams, prm::BCMLDTParams; kwargs...) =
    _solve_bcmldt_sue(min_cost_table_mt, bcmldt_loading_mt!, net, od, cost_params, prm; kwargs...)
