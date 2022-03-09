"""
`create_events`

Create an events matrix.

**Arguments**
- `events_description`: Events description as an array of tuples. For example,
    `[(:LL, 2*24), (:LD, 4, 12, 12), (:DD, 5*24)]` decodes 2 days (2 * 24
    hours) in constant light, 4 repetitions of a 12:12 light:dark cycle, and 5
    days (5 * 24 hourse) in constant darkness.

**Returns**
- `events`: Two-column matrix decoding times of light onsets (1st column) and
    light offsets (2nd columns).
"""
function create_events(events_description::Array)

    # Time counter
    current_time = 0.0

    # Initialize matrix of empty events
    events = Matrix{Float64}(undef, 0, 2)

    for x in events_description

        # Get code of the event
        code = Symbol(x[1])

        if code == :LD  # light-dark period

            # Iterate days
            for _ = 1:x[2]

                if !isempty(events) && events[end, 2] == current_time
                    # Extend an already existing ligt period
                    events[end, 2] = current_time + x[3]
                else
                    # Add a new light-dark cycle
                    events = [events; [current_time current_time + x[3]]]
                end

                # Increment the time counter by the duration of the added day
                current_time = current_time + x[3] + x[4]

            end

        elseif code == :DD  # constant darkness

            # Increment the time counter by the duration of the dark period
            duration = x[2]
            current_time += duration

        elseif code == :LL  # constant light

            if !isempty(events) && events[end, 2] == current_time
                # Extend an already existing ligt period
                events[end, 2] = current_time + x[2]
            else
                # Add a new light period
                events = [events; [current_time current_time + x[2]]]
            end

            # Increment the time counter the duration of the light period
            current_time += x[2]

        else  # unknown code
            msg = "Unknown events code `$code`! Should be `:LL`, `:LD`, or `:DD`."
            err = OscillatorPopulationError(msg)
            throw(err)
        end
    end

    return events

end


"""
`create_events`

Create an events matrix.

**Optional Arguments**
- `args...`: Event description. Examples for possible options:
    - `:LL, 24` for constant light of length 24
    - `:DD, 24` for constant darkness of length 24
    - `:LD, 4, 12, 12` for four repetitions of a 12:12 light:dark cycle

**Returns**
- `events`: Two-column matrix decoding times of light onsets (1st column) and
    light offsets (2nd columns).
"""
function create_events(args...)

    if isempty(args)
        events = create_events([])
    else
        events = create_events([args])
    end

    return events

end


"""
`create_events_cycle`

Create an events matrix for a regular light:dark cycle.

**Arguments**
- `duration`: Time duration to be covered by the cycle.
- `period`: Period of the cycle.

**Optional Arguments**
- `photoperiod`: Photoperiod of a cycle (fraction of the period spend in light).

**Returns**
- `events`: Two-column matrix decoding times of light onsets (1st column) and
    light offsets (2nd columns).
"""
function create_events_cycle(duration, period, photoperiod=0.5)

    if photoperiod == 1
        events = create_events(:LL, duration)    
    elseif photoperiod == 0
        events = create_events(:DD, duration)
    elseif 0 < photoperiod < 1
        light_duration = period*photoperiod
        night_duration = period*(1-photoperiod)
        n_repetitions = ceil(duration / period)
        events = create_events(:LD, n_repetitions, light_duration, night_duration)
    else
        msg = "Photoperiod is $photoperiod, but must be between 0 and 1!"
        err = OscillatorPopulationError(msg)
        throw(err)
    end
    return events

end


"""
`events_to_function(events::Matrix)`

Convert an events matrix to a function.

**Arguments**
- `events`: Events matrix.

**Returns**
- `fun`: A function that return 1 for light periods and 0 for dark periods.
"""
function events_to_function(events::Matrix)
    N = size(events, 1)
    function fun(t)
        for i = 1:N
            if events[i, 1] <= t < events[i, 2]
                return 1.0
            end
        end
        return 0.0
    end
    return fun
end


"""
`plot_events`

Plot events.

**Arguments**
- `events`: Events matrix.

**Keyword Arguments**
- `ax`: `PyPlot` axes.
- `ylims`: Range of y-axis values for the plotted events specified as `[y_min, y_max]`.
- `color`: Color of the plotted events.
- `zorder`: Order of the events in the `PyPlot` layers.

**Returns**
- `events`: Two-column matrix decoding times of light onsets (1st column) and
    light offsets (2nd columns).
"""
function plot_events(events; ax=gca(), ylims=[], color="#ffbfbfff", zorder=0)

    if isempty(ylims)
        # Get ylims based on the plotted lines of the Axes

        lines = ax.get_lines()
        if isempty(lines)
            # If there are no lines, default ylims = [0, 1]
            B = 0.0
            T = 1.0
        else
            # If there are lines, ylims = [min(lines), max(lines)]
            B = Inf
            T = -Inf
            for line in lines
                y = line.get_ydata()
                B = min(B, minimum(y))
                T = max(T, maximum(y))
            end
        end
    else
        B = ylims[1]
        T = ylims[2]
    end

    n = size(events, 1)
    for i = 1:n
        # Print a box for each event
        L = events[i, 1]
        R = events[i, 2]
        ax.fill_between([L, R], [B, B], [T, T], linewidth=0.0, color=color,
            zorder=zorder)
    end

    return ax

end


"""
`create_callback(events::Matrix, i)`

Generate `PresetTimeCallback` from `events` modyfing `i`-th parameter.

**Arguments**
- `events`: Events matrix.
- `i`: Index of the parameter to be modified by the callback.

**Keyword Arguments**
- `callback_type`: "continuous" or "discrete".
"""
function create_callback(events::Matrix, i; callback_type="continuous")

    p = NaN  # holds the original value of the parameter

    # Initialization at the beginning of the integration
    initialize = function (c, u, t, integrator)

        # Protect the original parameters from overwriting
        integrator.p = copy(integrator.p)

        # Save the default parameter
        p = integrator.p[i]

        # Set the parameter to 0
        integrator.p[i] = 0.0

    end

    if callback_type == "continuous"
    
        # Set the parameter to its default value at the beginning of the event
        tstops_on = events[:, 1]
        affect_on! = (integrator) -> integrator.p[i] = p
        callback_on = PresetTimeCallback(tstops_on, affect_on!,
            save_positions=(false, false), initialize=initialize)

        # Set the parameter to zero at the end of the event
        tstops_off = events[:, 2]
        affect_off! = (integrator) -> integrator.p[i] = 0.0
        callback_off = PresetTimeCallback(tstops_off, affect_off!,
            save_positions=(false, false))

        # Join the on and off callbacks
        callback = CallbackSet(callback_on, callback_off)
    
    elseif callback_type == "discrete"

        # Convert the events to a function
        fun = events_to_function(events)

        # Call the function at every step of the integration
        affect! = (integrator) -> integrator.p[i] = fun(integrator.t)*p
        condition = (u, t, integrator) -> true
        callback = DiscreteCallback(condition, affect!; initialize=initialize,
            save_positions=(false, false))

    else

        msg = "Unknown callback type `$callback_type`! Possible options are
            `continuous` and `discrete`."
        err = OscillatorPopulationError(msg)
        throw(err)

    end

    return callback

end
