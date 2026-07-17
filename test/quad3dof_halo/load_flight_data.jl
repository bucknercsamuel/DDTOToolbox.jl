"""
Load flight-test pickles directly with Pickle.jl (pure Julia, no PyCall / no
juliacall). All numeric data is decoded into Julia `Array{T,N}` with correct
shape and dtype by Pickle.jl's native numpy handlers.

Caveat: Pickle.jl ships with NumPy 1.x reducer names (`numpy.core.multiarray.*`).
NumPy 2.0 renamed those internals to `numpy._core.multiarray.*` (note the
leading underscore on `_core`), so pickles produced from NumPy 2.x leave every
array as a `Defer` wrapper. We register the new names below so the existing
`np_multiarray_reconstruct` / `np_scalar` decoders fire and we get raw arrays.

Usage:
    include("load_flight_data.jl")
    data = load_flight_pickle(path; inspect=true)   # prints type tree
    sim_state = data["sim_state"]                   # already a Julia Matrix
"""

using Pickle
using Statistics
using DataInterpolations

const Defer = Pickle.Defer  # mutable struct: .head::Symbol, .args::Vector{Any}


"""Build an `NpyPickler` and register NumPy 2.x reducer paths so post-NumPy-2.0
pickles decode their arrays/scalars natively instead of leaving them as
`Defer(:build, ...)` wrappers."""
function _flight_pickler()
    p = Pickle.NpyPickler()
    p.mt["numpy._core.multiarray._reconstruct"] = Pickle.np_multiarray_reconstruct
    p.mt["numpy._core.multiarray.scalar"]       = Pickle.np_scalar
    return p
end


# -----------------------------------------------------------------------------
# Inspection: print the type tree so we can see exactly what came in.
# -----------------------------------------------------------------------------
function _typestr(x; max_str_len = 60)
    if x isa AbstractArray
        return "$(typeof(x).name.name){$(eltype(x))}$(size(x))"
    elseif x isa Defer
        return "Defer(:$(x.head), $(length(x.args)) args)"
    elseif x isa AbstractDict
        return "Dict($(length(x)) keys)"
    elseif x isa AbstractString
        s = "\"$(x)\""
        return length(s) <= max_str_len ? s : "\"$(first(x, max_str_len-5))...\""
    elseif x isa Tuple
        return "Tuple($(length(x)) elems)"
    elseif x isa AbstractVector
        return "Vector{$(eltype(x))}($(length(x)))"
    else
        return string(typeof(x))
    end
end

"""Recursively print the structure of `x` (dicts, lists, Defers, arrays)."""
function inspect_pickle(x; max_depth = 4, _depth = 0, _label = "")
    pad = "  " ^ _depth
    label = isempty(_label) ? "" : "$_label : "
    if _depth >= max_depth
        println(pad, label, _typestr(x), "  [...max depth...]")
        return
    end
    if x isa AbstractDict
        println(pad, label, "Dict ($(length(x)) keys)")
        for (k, v) in x
            inspect_pickle(v; max_depth = max_depth, _depth = _depth + 1, _label = repr(k))
        end
    elseif x isa Defer
        println(pad, label, "Defer(:$(x.head), $(length(x.args)) args)")
        for (i, a) in enumerate(x.args)
            inspect_pickle(a; max_depth = max_depth, _depth = _depth + 1, _label = "arg[$i]")
        end
    elseif x isa AbstractArray && eltype(x) === Any
        n = length(x)
        println(pad, label, "Array{Any}$(size(x))  (object array, $n elems)")
        for i in 1:min(n, 3)
            inspect_pickle(x[i]; max_depth = max_depth, _depth = _depth + 1, _label = "[$i]")
        end
        if n > 3
            println(pad, "  ... ($(n - 3) more)")
        end
    elseif x isa Tuple
        println(pad, label, "Tuple($(length(x)) elems)")
        for (i, a) in enumerate(x)
            inspect_pickle(a; max_depth = max_depth, _depth = _depth + 1, _label = "[$i]")
        end
    elseif x isa AbstractVector && !(x isa AbstractVector{<:Number})
        n = length(x)
        println(pad, label, "Vector{$(eltype(x))}($n)")
        for i in 1:min(n, 3)
            inspect_pickle(x[i]; max_depth = max_depth, _depth = _depth + 1, _label = "[$i]")
        end
        if n > 3
            println(pad, "  ... ($(n - 3) more)")
        end
    else
        println(pad, label, _typestr(x))
    end
end


# -----------------------------------------------------------------------------
# Loading: pure Julia, no Python at all.
# -----------------------------------------------------------------------------
function load_flight_pickle(path; inspect = false, inspect_depth = 4)
    data = open(path, "r") do io
        Pickle.load(_flight_pickler(), io)
    end
    if inspect
        println("\n=== Structure of $(basename(path)) ===")
        inspect_pickle(data; max_depth = inspect_depth)
        println()
    end
    return data
end


# -----------------------------------------------------------------------------
# Post-processing: remove time discontinuities and zero the start.
# -----------------------------------------------------------------------------
"""
Slice a flight `data` dict by a manually-specified list of bad sample indices,
then enforce time monotonicity, then zero the time start.

Pipeline:

1. Drop any indices listed in `exclude_indices` (1-based into the original
   `sim_time`).
2. Walk the surviving indices left-to-right and enforce `sim_time` is
   non-decreasing by removing any index that breaks it: whenever the most
   recently kept index has a strictly larger `sim_time` than the current
   candidate, pop the kept index. This implements "if `time[idx] >
   time[idx+1]`, remove `idx`" iteratively until no violations remain.
3. Slice every per-sample array (matrix with `size(·, 2) == n_orig`, vector
   with `length(·) == n_orig`) by the surviving original indices.
4. Subtract `sim_time[first_kept]` from `sim_time` and `guid_update_times` so
   the cleaned trace starts at `t = 0`.

Returns a *new* dict (input is not mutated).
"""
function postprocess_flight_data(data; exclude_indices = Int[])
    sim_time_orig = vec(collect(data["sim_time"]))
    n_orig = length(sim_time_orig)

    exclude_set = Set(exclude_indices)
    out_of_range = filter(i -> i < 1 || i > n_orig, collect(exclude_set))
    if !isempty(out_of_range)
        @warn "postprocess_flight_data: exclude_indices contains out-of-range values; ignoring them" out_of_range n_orig
    end
    candidate_indices = filter(i -> !(i in exclude_set), 1:n_orig)

    # Enforce non-decreasing sim_time. For each new candidate, pop any
    # previously-kept indices whose sim_time exceeds the current one.
    kept_indices = Int[]
    unsorted_removed = Int[]
    for i in candidate_indices
        while !isempty(kept_indices) && sim_time_orig[kept_indices[end]] > sim_time_orig[i]
            push!(unsorted_removed, pop!(kept_indices))
        end
        push!(kept_indices, i)
    end

    if isempty(kept_indices)
        @warn "postprocess_flight_data: every sample was excluded; returning unchanged"
        return data
    end

    new_data = Dict{Any, Any}()
    for (key, val) in data
        if val isa AbstractMatrix && size(val, 2) == n_orig
            new_data[key] = val[:, kept_indices]
        elseif val isa AbstractVector && length(val) == n_orig
            new_data[key] = val[kept_indices]
        else
            new_data[key] = val
        end
    end

    t_shift = new_data["sim_time"][1]
    new_data["sim_time"] = new_data["sim_time"] .- t_shift
    if haskey(new_data, "guid_update_times")
        new_data["guid_update_times"] = vec(collect(new_data["guid_update_times"])) .- t_shift
    end

    sort!(unsorted_removed)
    preview_unsorted = unsorted_removed[1:min(10, length(unsorted_removed))]
    @info "postprocess_flight_data" n_total=n_orig n_kept=length(kept_indices) n_excluded_manual=length(exclude_set) n_unsorted_removed=length(unsorted_removed) excluded=sort(collect(exclude_set)) first_unsorted_removed=preview_unsorted t_shift=round(t_shift; digits = 4)
    return new_data
end


"""
Close large time gaps in a flight `data` dict by shifting time forward.

For every pair of consecutive samples where `sim_time[i+1] - sim_time[i] >
threshold`, the full gap (the `dt` value at that step) is subtracted from
`sim_time` at all indices `>= i+1`. The samples themselves are kept — only
their timestamps move. After processing, the offending `dt` becomes 0 (the two
samples share a timestamp) and every subsequent timestamp is shifted earlier
by the same amount.

`guid_update_times` is adjusted in parallel: each guidance update time has
the total duration of all gaps that fell strictly before it (in the original
`sim_time`) subtracted, so it stays aligned with the data.

Returns a new dict (input is not mutated). An `@info` line reports how many
gaps were closed, the total time removed, and where the gaps were.
"""
function compress_time_gaps(data; threshold = 0.1)
    sim_time_orig = vec(collect(data["sim_time"]))
    n = length(sim_time_orig)
    if n < 2
        return data
    end

    new_sim_time = copy(sim_time_orig)
    gap_records = Tuple{Int, Float64}[]  # (i, dt_removed) for each closed gap

    for i in 1:n-1
        dt = new_sim_time[i+1] - new_sim_time[i]
        if dt > threshold
            new_sim_time[i+1:end] .-= dt
            push!(gap_records, (i, dt))
        end
    end

    new_data = Dict{Any, Any}(k => v for (k, v) in data)
    new_data["sim_time"] = new_sim_time

    if haskey(new_data, "guid_update_times")
        guid_times_orig = vec(collect(new_data["guid_update_times"]))
        adjusted = similar(guid_times_orig)
        for k in eachindex(guid_times_orig)
            shift = 0.0
            for (i, dt_removed) in gap_records
                # Gap at step i was originally between sim_time_orig[i] and sim_time_orig[i+1].
                # If guidance update happened at-or-after the post-gap sample, subtract.
                if guid_times_orig[k] >= sim_time_orig[i+1]
                    shift += dt_removed
                end
            end
            adjusted[k] = guid_times_orig[k] - shift
        end
        new_data["guid_update_times"] = adjusted
    end

    total_removed = isempty(gap_records) ? 0.0 : sum(last, gap_records)
    preview = gap_records[1:min(5, length(gap_records))]
    @info "compress_time_gaps" n=n n_gaps_closed=length(gap_records) total_time_removed=round(total_removed; digits = 4) threshold=threshold first_gaps=preview
    return new_data
end


"""
Resample a flight `data` dict onto a uniform time grid with step `dt` (seconds).
`sim_time`, `sim_state`, and `sim_control` (if present and same length as
`sim_time`) are resampled column-by-column; all other fields pass through
unchanged.

`interpolation` selects the interpolation kind:

* `:cubic`  — `DataInterpolations.CubicSpline` (smooth, can overshoot across
              large gaps)
* `:linear` — `DataInterpolations.LinearInterpolation` (no overshoot, but
              kinked at every sample)
* `:hybrid` — `alpha * cubic + (1 - alpha) * linear`, evaluated pointwise. With
              `alpha = 1` this matches `:cubic`; with `alpha = 0` it matches
              `:linear`. Use intermediate values (e.g. `alpha = 0.5`) to damp
              cubic-spline overshoot near large gaps while preserving smoothness
              away from them.

`alpha` is ignored unless `interpolation = :hybrid`. It must lie in `[0, 1]`.

If `sim_time` contains duplicate timestamps (e.g. from `compress_time_gaps`),
they are nudged apart by `1e-9 s` so the interpolator sees strictly increasing
time. If `sim_time` is not sorted, the function aborts with a warning — apply
`postprocess_flight_data` first.
"""
function resample_flight_data(data; dt = 0.05, interpolation::Symbol = :cubic, alpha::Real = 0.5)
    if !(interpolation in (:cubic, :linear, :hybrid))
        throw(ArgumentError("resample_flight_data: interpolation must be :cubic, :linear, or :hybrid; got :$interpolation"))
    end
    if interpolation === :hybrid && !(0 <= alpha <= 1)
        throw(ArgumentError("resample_flight_data: alpha must be in [0, 1], got $alpha"))
    end

    sim_time = vec(collect(data["sim_time"]))
    n = length(sim_time)
    if n < 4
        @warn "resample_flight_data: fewer than 4 samples; returning unchanged"
        return data
    end
    if !issorted(sim_time)
        @warn "resample_flight_data: sim_time is not sorted; returning unchanged. Apply postprocess_flight_data first."
        return data
    end

    t_fit = copy(sim_time)
    n_nudged = 0
    for i in 2:n
        if t_fit[i] <= t_fit[i-1]
            t_fit[i] = t_fit[i-1] + 1e-9
            n_nudged += 1
        end
    end

    t_grid = collect(range(sim_time[1], sim_time[end]; step = dt))

    function interp_rows(matrix)
        out = Array{Float64}(undef, size(matrix, 1), length(t_grid))
        for k in 1:size(matrix, 1)
            y = Vector{Float64}(matrix[k, :])
            if interpolation === :cubic
                itp = CubicSpline(y, t_fit)
                out[k, :] = itp.(t_grid)
            elseif interpolation === :linear
                itp = LinearInterpolation(y, t_fit)
                out[k, :] = itp.(t_grid)
            else  # :hybrid
                itp_c = CubicSpline(y, t_fit)
                itp_l = LinearInterpolation(y, t_fit)
                @inbounds for j in eachindex(t_grid)
                    out[k, j] = alpha * itp_c(t_grid[j]) + (1 - alpha) * itp_l(t_grid[j])
                end
            end
        end
        return out
    end

    new_data = Dict{Any, Any}(k => v for (k, v) in data)
    new_data["sim_time"] = t_grid

    if haskey(new_data, "sim_state")
        m = new_data["sim_state"]
        if m isa AbstractMatrix && size(m, 2) == n
            new_data["sim_state"] = interp_rows(m)
        else
            @warn "resample_flight_data: sim_state size doesn't match sim_time; not resampled" size=size(m) n=n
        end
    end

    if haskey(new_data, "sim_control")
        m = new_data["sim_control"]
        if m isa AbstractMatrix && size(m, 2) == n
            new_data["sim_control"] = interp_rows(m)
        else
            @warn "resample_flight_data: sim_control size doesn't match sim_time; not resampled" size=size(m) n=n
        end
    end

    @info "resample_flight_data" n_orig=n n_new=length(t_grid) dt=dt interpolation=interpolation alpha=(interpolation === :hybrid ? alpha : nothing) time_span=(round(sim_time[1]; digits = 4), round(sim_time[end]; digits = 4)) duplicates_nudged=n_nudged
    return new_data
end
