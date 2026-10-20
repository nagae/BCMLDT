# =====================================================================
#  Reading and downloading TNTP files (https://github.com/bstabler/TransportationNetworks)
# =====================================================================

# If `ln` is a metadata line `<KEY> value`, store it in `meta` and return true.
function _tntp_metadata!(meta::Dict{String,String}, ln::AbstractString)
    m = match(r"^\s*<(.+?)>\s*(.*)$", ln)
    m === nothing && return false
    meta[uppercase(strip(m.captures[1]))] = strip(m.captures[2])
    return true
end

"""
    read_tntp_net(path) -> (net, cost_params, len)

Read a `*_net.tntp` file. Returns the `Network`, the `BPRParams` with `t_add = 0`,
and the link lengths `len` (in the unit of the file), all in the order of the links
in the file.
"""
function read_tntp_net(path::AbstractString)
    meta = Dict{String,String}()
    rows = Vector{Vector{Float64}}()
    body = false
    for ln in eachline(path)
        if !body
            _tntp_metadata!(meta, ln) && continue
            if startswith(ln, "~")
                body = true
            end
            continue
        end
        str = strip(ln)
        (isempty(str) || startswith(str, "~")) && continue
        v = split(replace(str, ";" => " "))
        push!(rows, parse.(Float64, v[1:7]))       # tail, head, capacity, length, free flow time, b, power
    end
    tbl = reduce(vcat, permutedims.(rows))        # one row per link
    num_nodes = parse(Int, meta["NUMBER OF NODES"])
    ftn = parse(Int, get(meta, "FIRST THRU NODE", "1"))
    net = Network(num_nodes, Int.(tbl[:, 1]), Int.(tbl[:, 2]); first_thru_node = ftn)
    cost_params = BPRParams(tbl[:, 5], tbl[:, 3], tbl[:, 6], tbl[:, 7], zeros(size(tbl, 1)))   # free flow time, capacity, b, power
    return (net = net, cost_params = cost_params, len = tbl[:, 4])
end

"""
    read_tntp_trips(path) -> ODPairs

Read a `*_trips.tntp` file. OD pairs with zero demand and with `r == s` are dropped,
and the pairs are ordered by origin.
"""
function read_tntp_trips(path::AbstractString)
    txt = read(path, String)
    meta = Dict{String,String}()
    head, body = split(txt, "<END OF METADATA>"; limit = 2)
    foreach(ln -> _tntp_metadata!(meta, ln), split(head, '\n'))
    nz = parse(Int, meta["NUMBER OF ZONES"])
    Q = zeros(nz, nz)
    for blk in split(body, r"Origin\s+")[2:end]
        o_str, rest = split(blk, '\n'; limit = 2)
        o = parse(Int, strip(o_str))
        for m in eachmatch(r"(\d+)\s*:\s*([\d.eE+-]+)", rest)
            Q[o, parse(Int, m.captures[1])] = parse(Float64, m.captures[2])
        end
    end
    return ODPairs(Q)
end

"""
    read_tntp_ue(path, net) -> (flow, cost)

Read a `*_flow.tntp` file (a UE solution provided with the data) and return the
link flows and link costs in the order of the links of `net`. Note that this file
is not provided for every network in the repository.
"""
function read_tntp_ue(path::AbstractString, net::Network)
    isfile(path) || error("$(path) not found: a UE solution is not provided for every network")
    id = Dict{Tuple{Int,Int},Int}()
    for a in 1:num_links(net)
        key = (net.tail[a], net.head[a])
        haskey(id, key) && error("parallel links $(key) cannot be matched to the flow file")
        id[key] = a
    end
    flow = fill(NaN, num_links(net))
    cost = fill(NaN, num_links(net))
    for ln in Iterators.drop(eachline(path), 1)    # header: From To Volume Cost
        v = split(ln)
        length(v) < 4 && continue
        a = id[(parse(Int, v[1]), parse(Int, v[2]))]
        flow[a] = parse(Float64, v[3])
        cost[a] = parse(Float64, v[4])
    end
    any(isnan, cost) && error("some links are missing in $(path)")
    return (flow = flow, cost = cost)
end

# ---------------------------------------------------------------------
#  Downloading TNTP files
# ---------------------------------------------------------------------

# Commit of the TransportationNetworks repository used in the paper (2 August 2023)
const TNTP_COMMIT = "375e0da93858c547230c5cf9ea8a96de4ccff29e"
# Folders in the repository whose names differ from the file names
const _TNTP_FOLDER = Dict("ChicagoSketch" => "Chicago-Sketch")

"""
    download_tntp(name, dir) -> dir

Download `name_net.tntp`, `name_trips.tntp`, and `name_flow.tntp` (if provided) of the
network `name` (e.g. `"Anaheim"`) from the TransportationNetworks repository
(https://github.com/bstabler/TransportationNetworks) at the commit used in the paper,
and save them in `dir`. Files already in `dir` are not downloaded again.
"""
function download_tntp(name::AbstractString, dir::AbstractString)
    mkpath(dir)
    folder = get(_TNTP_FOLDER, name, name)
    base = "https://raw.githubusercontent.com/bstabler/TransportationNetworks/$(TNTP_COMMIT)/$(folder)"
    for kind in ("net", "trips", "flow")
        file = joinpath(dir, "$(name)_$(kind).tntp")
        isfile(file) && continue
        try
            Downloads.download("$(base)/$(name)_$(kind).tntp", file)
        catch err
            kind == "flow" || rethrow(err)                 # a UE solution is not provided for every network
        end
    end
    return dir
end

