"""
`rsquared`

Calculate the coefficient of determination.

**Arguments**
- `x`: Vector of real data.
- `y`: Vector of predicted values.

**Returns**
- `R`: The coefficient of determination as a number on the interval `[-∞, 1]`.

**Reference**
- https://doi.org/10.1016/0022-1694(70)90255-6
"""
function rsquared(x, y)

    μ = mean(x)
    SSres = sum((x .- y).^2)
    SStot = sum((x .- μ).^2)
    R = 1 - SSres/SStot
    return R

end


"""
`cmean`

Compute circular mean of an array.

**Arguments**
- `x`: Array of angles in the interval [0, 1].

**Returns**
- `m`: Circular mean of `x`.
"""
function cmean(x)
    x = mod.(x, 1)
    x_rad = x * 2π
    S = mean(sin.(x_rad))
    C = mean(cos.(x_rad))
    m_rad = atan(S, C) 
    m = m_rad / 2π
    m = mod(m, 1)
    return m
end


"""
`cstd(x)`

Compute circular standard deviation of an array.

**Arguments**
- `x`: Array of angles normalized to the interval [0, 1].

**Returns**
- `sd`: Circular standard deviation of `x`.
"""
function cstd(x)
    x = mod.(x, 1)
    x_rad = x * 2π
    S = mean(sin.(x_rad))
    C = mean(cos.(x_rad))
    R = sqrt(S^2 + C^2)
    x = -2*log(R)
    if x < 0  # for numerical errors close to 0
        x = 0
    end
    sd_rad = sqrt(x)
    sd = sd_rad / 2π
    return sd
end


"""
`window_xcorr`

Compute moving-window cross-correlation. This is a particular implementation of
the cross-correlation function that assumes that `x` is shorter than `y`. `x` is
moved along `y` and the cross-correlation function is calculated as weighted
moving average, where `x` are the weights and `y` is the averaged signal.

**Arguments**
- `x`: A signal array.
- `y`: A signal array.

**Returns**
- `r`: Cross-correlation of `x` and `y`.
"""
function window_xcorr(x, y)

    # Remove mean
    x = zscore(x)
    y = zscore(y)

    # Length of the input and output vectors
    nx = length(x)
    ny = length(y)
    nr = ny - nx + 1

    # Check that x is shorter than y
    if nx > ny
        msg = "x must be shorter than y!"
        err = OscillatorPopulationError(msg)
        throw(err)
    end

    xx = vcat(x, fill(0, ny-nx))
    r = crosscov(xx, y, 0:(nr-1); demean=false) .* ny ./ nx

    return r

end


"""
`cxcorr`

Compute circular cross-correlation.

**Arguments**
- `x`: Array representing one period of a signal.
- `y`: Array representing one period of a signal.

**Returns**
- `r`: Circular cross-correlation of `x` and `y`.
"""
function cxcorr(x, y)

    x = zscore(x)
    y = zscore(y)

    nx = length(x)
    xx = vcat(x, fill(0, nx))
    yy = vcat(y, y)
    r = crosscov(xx, yy, 0:(nx-1); demean=false) .* length(xx) ./ nx

    return r

end


"""
`estimate_phase_array`

Estimate entrainment phase. The input parameters are expected to represent a
model output entrained to a periodic square input. The first and last event are
discarted.

**Arguments**
- `t`: Time vector.
- `x`: Data vector or matrix, where each column represents one data vector.
- `events`: Events matrix with two columns representing a square signal.

**Keyword Arguments**
- `method`: Choose method to estimate the phase. Possible options are 
    `"peak_height"` (default), `"peak_prominence"`, and `"cxcorr"`.
- `smooth_span`: Length of the moving average filter (in samples) used to average
    `x` before the phase is estimated. Default value is 1 (no smoothing).
- `show_plots`: [NOT IMPLEMENTED]

**Returns**
- `phase_array`: Array of phases estimated at each cycle of the input data. If 
    `x` is a matrix, `phase_array` is also a matrix, where each column contains
    estimated phases for the individual data vectors.
"""
function estimate_phase_array(t, x, events; method="peak_height",
    smooth_span=1, show_plots=false)

    if smooth_span > 1
        x = smooth(x, span=smooth_span)
    end

    if method == "peak_height"
        phase_array = estimate_phase_array_peaks(t, x, events; use_prominence=false)
    elseif method == "peak_prominence"
        phase_array = estimate_phase_array_peaks(t, x, events; use_prominence=true)
    elseif method == "cxcorr"
        phase_array = estimate_phase_array_cxcorr(t, x, events)
    end

    return phase_array

end


"""
`estimate_phase_array_peaks`

Estimate entrainment phase based on the peak prominence.

**Arguments**
- `t`: Time vector.
- `x`: Data vector or matrix, where each column represents one data vector.
- `events`: Events matrix with two columns representing a square signal.

**Keyword Arguments**
- `use_prominence`: If `true`, the algorithm will use peak prominences
    instead of peak heights to determine the heights peaks. Default value is
    `false`.

**Returns**
- `phase_array`: Array of phases estimated at each cycle of the input data.
"""
function estimate_phase_array_peaks(t, x, events; use_prominence=false)

    pr = findpeaks(x, t)
    if use_prominence
        pks = peakprominences(pr)
    else
        pks = peakheights(pr)
    end
    locs = peaklocations(pr)

    n_events = size(events, 1)
    phase_array = Float64[]
    for i_event = 2:(n_events-1)
        window_start = events[i_event, 1]
        window_end = events[i_event+1, 1]
        idx = window_start .<= locs .< window_end
        window_locs = locs[idx]
        window_pks = pks[idx]

        if isempty(window_pks)
            # If there are no peaks in the window, set phase to NaN
            push!(phase_array, NaN)
        else
            # If there are peaks in the window, estimate phase
            idx_peak = argmax(window_pks)
            loc = window_locs[idx_peak]
            phase = loc - window_start
            phase /= (window_end - window_start)
            push!(phase_array, phase)
        end
    end

    return phase_array

end


"""
`estimate_phase_array_cxcorr`

Estimate entrainment phase based on the circular cross-correlation.

**Arguments**
- `t`: Time vector.
- `x`: Data vector or matrix, where each column represents one data vector.
- `events`: Events matrix with two columns representing a square signal.

**Returns**
- `phase_array`: Array of phases estimated at each cycle of the input data.
"""
function estimate_phase_array_cxcorr(t, x, events)

    fun = events_to_function(events)
    y = fun.(t)

    n_events = size(events, 1)
    phase_array = Float64[]
    for i_event = 2:(n_events-1)
        window_start = events[i_event, 1]
        window_end = events[i_event+1, 1]
        window_idices = window_start .<= t .< window_end
        if sum(window_idices) == 0
            push!(phase_array, NaN)
        else
            r = cxcorr(x[window_idices], y[window_idices])
            _, max_index = findmax(r)
            phase = (max_index - 1) / sum(window_idices)
            push!(phase_array, phase)
        end
    end

    return phase_array

end


function estimate_phase_array(t, X::Matrix, events; kwargs...)

    # Call estimate_phase_array on each column of a matrix
    n_trajectories = size(X, 2)
    n_events = size(events, 1) - 2
    phase_matrix = fill(NaN, n_events, n_trajectories)
    for i_trajectory = 1:n_trajectories
        x = X[:, i_trajectory]
        phase_array = estimate_phase_array(t, x, events; kwargs...)
        phase_matrix[:, i_trajectory] = phase_array
    end
    return phase_matrix

end


"""
`estimate_order_parameter`

Estimate period of a signal based on its autocorrelation function.

**Argument**
- `phase_array`: Array of phases normalized to interval [0, 1].

**Keyword Arguments**
- `skip_nan`: Skip NaN values. If `false`, `NaN` values are propageted to the
    output. Default value is `true`.

**Returns**
- `phase_coherence`: A number between 0 (no coherence) to 1 (full coherence).
    If `phase_array` is a matrix, phase coherence is calculated as a mean of
    the phase coherences calculated from individual rows.
- `collective_phase`: A number between 0 and 1 indicating the average phase.
    If `phase_array` is a matrix, collective phase is circular mean of the
    indivudal phases.
"""
function estimate_order_parameter(phase_array; skip_nan=true)

    if minimum(phase_array) < 0 || maximum(phase_array) > 1
        msg = "Phase must be on interval [0, 1]!"
        err = OscillatorPopulationError(msg)
        throw(err)
    end

    nan_indices = isnan.(phase_array)
    if !skip_nan && any(nan_indices)
        return NaN, NaN
    else
        phase_array = phase_array[.!nan_indices]
        if isempty(phase_array)
            return NaN, NaN
        end
    end

    n = length(phase_array)
    complex_order_parameter = 0im
    for phase in phase_array
        complex_order_parameter += exp(2π * phase * im)
    end
    complex_order_parameter /= n
    phase_coherence = abs(complex_order_parameter)
    collective_phase = angle(complex_order_parameter)  # ∈ [-π, π]
    if collective_phase < 0  # map from [-π, 0] on [π, 2π]
        collective_phase += 2π
    end
    collective_phase /= 2π  # map from [0, 2π] on [0, 1]
    return phase_coherence, collective_phase

end


function estimate_order_parameter(phase_matrix::Matrix)

    # Call estimate_order_parameter on each row of a matrix
    n = size(phase_matrix, 1)
    phase_coherence_arr = fill(NaN, n)
    collective_phase_arr = fill(NaN, n)
    for i = 1:n
        phase_array = phase_matrix[i, :]
        phase_coherence, collective_phase = estimate_order_parameter(phase_array)
        phase_coherence_arr[i] = phase_coherence
        collective_phase_arr[i] = collective_phase
    end
    phase_coherence_mean = mean(phase_coherence_arr)
    collective_phase_mean = cmean(collective_phase_arr)
    return phase_coherence_mean, collective_phase_mean

end


"""
`estimate_period(t, x; kwargs...)`

Estimate period of a signal based on its autocorrelation function.

**Argument**
- `t`: Time vector.
- `x`: Data vector.

**Keyword Arguments**
- `n_lags`: Number of lags used for autocorrelation. Default is `length(x)`.
- `show_plots`: Visualize the result. Default is false.

**Returns**
- `period`: Estimated period (the position of the first peak).
- `peak`: Height of the first peak.
"""
function estimate_period(t, x; n_lags=length(x), show_plots=false)

    # If x is constat, return NaNs
    if all(x[1] .== x)
        return (NaN, NaN)
    end

    # Normalize input
    x = zscore(x)

    # Calculate autocorrelation
    n = min(n_lags, length(x))  # number of lags to calculate
    R = autocor(x, 0:(n-1))  # autocorrelation function

    # Find peaks of the autocorrelation function
    pr = findpeaks(R, (0:(n-1)) .* mean(diff(t)))  # , sortstr="descend", sortref="prominence")

    # If there are no peaks, return NaNs
    if length(pr) == 0
        return (NaN, NaN)
    end

    # Get the first peak and its location as the estimated period
    peak = peakheights(pr)[1]
    period = peaklocations(pr)[1]

    # Visualize the results
    if show_plots
        fig, axs = subplots(2)
        axs[1].plot(t, x, color="black")
        axs[1].set_title("Time Series")
        axs[2].plot((0:(n-1)) .* mean(diff(t)), R, color="black")
        axs[2].plot([period, period], [0, peak], color="red")
        axs[2].set_title("peak = $(round(peak, digits=2)), " *
            "period = $(round(period, digits=2))")
        fig.tight_layout()
    end

    return (period, peak)

end


"""
`estimate_winding_number`

Estimate the winding number around the origin.

**Argument**
- `x`: Data vector for the x-coordinate.
- `y`: Data vector for the y-coordinate.

**Keyword Argument**
- `remove_mean`: Subtract mean from `x` and `y` before estimating the winding
    number. Default value is `true`.

**Returns**
- `winding_number`: Estimated winding number.
"""
function estimate_winding_number(x, y; remove_mean=true)

    if remove_mean
        x = x .- mean(x)
        y = y .- mean(y)
    end
    n = length(x)
    total_angle = 0
    for i = 1:n-1
        total_angle += atan(x[i]*y[i+1] - y[i]*x[i+1], x[i]*x[i+1] + y[i]*y[i+1])
    end
    winding_number = abs(total_angle / 2π)
    return winding_number

end


"""
`estimate_period_winding_number`

Estimate period based on the winding number.

**Argument**
- `x`: Data vector for the x-coordinate.
- `y`: Data vector for the y-coordinate.
- `time_duration`: Time duration of `x` and `y`.

**Keyword Argument**
- `remove_mean`: Subtract mean from `x` and `y` before estimating the winding
    number. Default value is `true`.

**Returns**
- `period`: Estimated period.
"""
function estimate_period_winding_number(x, y, time_duration; remove_mean=true)

    winding_number = estimate_winding_number(x, y; remove_mean=remove_mean)
    period = time_duration / winding_number
    return period

end
