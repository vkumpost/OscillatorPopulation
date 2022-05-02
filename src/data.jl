const BIOLUM_LOCAL_PATH = joinpath(@__DIR__,"..", "data", "biolumdata.xlsx")


"""
`load_biolum`

Load data from a bioluminescence experiment.

**Argument**
- `plate`: Name of the plate to load. Available plates are: `Plate U1`, 
    `Plate U2`, `Plate U3`, `Plate U4`, `Plate D1 A`, `Plate D2 A`,
    `Plate D1 B`, `Plate D2 B`.

**Keyword Arguments**
- `first_event`: The first event to be considered (default `0`).
- `last_event`: The last event to be considered (default `Inf`).
- `max_time`: Maximal time from the beginning or the first select event (if
    `first_event` was passed) to be loaded (default `Inf`).

**Returns**
- `df`: `DataFrame` with the loaded data. The first column `Time` contains
    timestamps, the second column `Light` is a binary vector with `1` indicating
    lights on and `0` lights off. The remaining columns are the recorded values
    of bioluminescence.
"""
function load_biolum(plate; first_event=0, last_event=Inf, max_time=Inf)

    table, header = XLSX.readtable(BIOLUM_LOCAL_PATH, plate)
    float_array = convert(Array{Array{Float64, 1}, 1}, table)
    df = DataFrame(float_array, header)

    # If last_event is passed, remove tailing samples
    if last_event < Inf

        events = _detect_events(df)
        t_end = events[last_event, 1]
        idx = df[:, "Time"] .< t_end
        df = df[idx, :]

    end

    # If first_event is passed, remove heading samples
    if first_event > 0

        events = _detect_events(df)
        t_start = events[first_event, 1]
        idx = t_start .< df[:, "Time"]
        df = df[idx, :]
        df[:, "Time"] .-= t_start

    end

    # If max_time is passed, remove tailing samples
    if max_time < Inf

        idx = df[:, "Time"] .<= max_time
        df = df[idx, :]

    end

    return df

end


"""
`_detect_events`

Detect events from a DataFrame.

**Arguments**
- `t`: Vector of timestamps.
- `x`: A binary vector.

**Returns**
- `events`: Events matrix representing event starts (first column) and event
    ends (second column).
"""
function _detect_events(t, x)

    # If there is no light recorded, return an empty matrix
    if all(x .== 0)
        return Matrix{Float64}(undef, 0, 2)
    end

    # Iterate vector x and find all event starts (x go from 0 to 1) and all
    # event ends (x go from 1 to 0).
    event_starts = Array{Float64, 1}()
    event_ends = Array{Float64, 1}()
    for i = 1:length(x)-1
        if x[i] < 0.5 && x[i+1] > 0.5
            push!(event_starts, (t[i] + t[i+1])/2)
        elseif x[i] > 0.5 && x[i+1] < 0.5
            push!(event_ends, (t[i] + t[i+1])/2)
        end
    end

    # Correct special cases for the first and last event
    if isempty(event_starts)
        event_starts = [t[1]; event_starts]
    end
    if isempty(event_ends)
        event_ends = [event_ends; t[end]]
    end
    if event_ends[1] < event_starts[1]
        event_starts = [t[1]; event_starts]
    end
    if event_ends[end] < event_starts[end]
        event_ends = [event_ends; t[end]]
    end

    # Build and return the events matrix
    events = [event_starts event_ends]
    return events
    
end


"""
`_detect_events`

Detect events from a DataFrame.

**Arguments**
- `df`: `DataFrame` with columns `Time` and `Light`.

**Returns**
- `events`: Events matrix representing event starts (first column) and event
    ends (second column).
"""
function _detect_events(df::DataFrame)

    # Extract time and light from DataFrame
    t = df[:, "Time"]
    x = df[:, "Light"]

    # Find and returns events
    events = _detect_events(t, x)
    return events
    
end


"""
`biolum_zscore_traces`

Apply zscore to the individual traces of bioluminescence data.

**Arguments**
- `df`: `DataFrame` with the bioluminescence data. First two columns are `Time`
    and `Light`.

**Returns**
- `df`: The copy of the input dataframe with zscored trajectories.

Apply function `fun` to each column of DataFrame `df` except for columns `Time`
and `Light`.
"""
function biolum_zscore_traces(df::DataFrame)
    df = deepcopy(df)
    for name in names(df)
        if name ∉ ("Time", "Light")
            df[:, name] = zscore(df[:, name])
        end
    end
    return df
end
