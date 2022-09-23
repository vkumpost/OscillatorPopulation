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
- `first_event_end`: If `true`, use the end instead of the beginning of the
    `first_event` to align the data (default `false`).
- `last_event_end`: If `true`, use the end instead of the beginning of the
    `end_event` to align the data (default `false`).
- `max_time`: Maximal time from the beginning or the first select event (if
    `first_event` was passed) to be loaded (default `Inf`).

**Returns**
- `df`: `DataFrame` with the loaded data. The first column `Time` contains
    timestamps, the second column `Light` is a binary vector with `1` indicating
    lights on and `0` lights off. The remaining columns are the recorded values
    of bioluminescence.
"""
function load_biolum(plate; first_event=0, last_event=Inf,
    first_event_end=false, last_event_end=false, max_time=Inf)
    # rename max_time to duration

    table, header = XLSX.readtable(BIOLUM_LOCAL_PATH, plate)
    float_array = convert(Array{Array{Float64, 1}, 1}, table)
    df = DataFrame(float_array, header)

    # If last_event is passed, remove tailing samples
    if last_event < Inf

        events = detect_events(df)
        if last_event_end
            t_end = events[last_event, 2]
        else
            t_end = events[last_event, 1]
        end
        idx = df[:, "Time"] .< t_end
        df = df[idx, :]

    end

    # If first_event is passed, remove heading samples
    if first_event > 0

        events = detect_events(df)
        if first_event_end
            t_start = events[first_event, 2]
        else
            t_start = events[first_event, 1]
        end
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
`detect_events`

Detect events from a DataFrame.

**Arguments**
- `df`: `DataFrame` with columns `Time` and `Light`.

**Returns**
- `events`: Events matrix representing event starts (first column) and event
    ends (second column).
"""
function detect_events(df::DataFrame)

    # Extract time and light from DataFrame
    t = df[:, "Time"]
    x = df[:, "Light"]

    # Find and returns events
    events = detect_events(t, x)
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


"""
`biolum_mean`

Calculate mean and standard deviation of the bioluminescence traces.

**Arguments**
- `df`: `DataFrame` with the bioluminescence data. First two columns are `Time`
    and `Light`.

**Returns**
- `df_mean`: `DataFrame` with columns `Time`, `Light`, `Mean`, and `STD`.
"""
function biolum_mean(df::DataFrame)
    df = deepcopy(df)
    df_mean = DataFrame(
        Time = df[:, "Time"],
        Light = df[:, "Light"],
        Mean = mean(Matrix(df[:, 3:end]), dims=2)[:, 1],
        STD = std(Matrix(df[:, 3:end]), dims=2)[:, 1]
    )
    return df_mean
end


"""
`save_data`

Save a dataframe as a csv file.

**Arguments**
- `df`: `DataFrame`.
- `filename`: Path to the csv file.

**Keyword Arguments**
- `force`: If `true`, `filename` will be overwritten, if already exists.
"""
function save_data(df::DataFrame, filename; force=false)

    # Check if the file already exists
    if !force && isfile(filename)
        msg = "File $(filename) already exists!"
        err = OscillatorPopulationError(msg)
        throw(err)
    end

    # Create necessary directories
    dir, _ = splitdir(filename)
    if !isdir(dir)
        mkpath(dir)
    end

    # Write the dataframe into the csv file
    matrix = Matrix(df)
    head = names(df)
    csv_content = [reshape(head, 1, length(head)); matrix]
    open(filename, "w") do io
        writedlm(io, csv_content, ',')
    end

end


"""
`load_data`

Load a dataframe from a csv file.

**Arguments**
- `filename`: Path to the csv file.
"""
function load_data(filename)

    matrix, head = readdlm(filename, ',', header=true)
    df = DataFrame(matrix, vec(head))
    return df

end


"""
`save_figure`

Save figure into a file.

**Arguments**
- `fig`: PyPlot figure.
- `filename`: Path to file.

**Keyword Arguments**
- `force`: If `true`, `filename` will be overwritten, if already exists.
"""
function save_figure(fig, filename; force=false)

    # Check if the file already exists
    if !force && isfile(filename)
        msg = "File $(filename) already exists!"
        err = OscillatorPopulationError(msg)
        throw(err)
    end

    # Create necessary directories
    dir, _ = splitdir(filename)
    if !isdir(dir)
        mkpath(dir)
    end

    # Write the dataframe into the csv file
    fig.savefig(filename)

end
