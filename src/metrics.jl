"""
`estimate_phase`

Estimate entrainment phase. The input parameters are expected to represent a
model output entrained to a periodic square input. The first and last event are
discarted.

**Arguments**
- `t`: Time vector.
- `x`: Data vector.
- `events`: Events matrix with two columns representing a square signal.

**Keyword Arguments**
- `normalize`: If `true`, the estimated phase is normalized by the event period.
- `show_plots`: If `true`, a figure is generated for a visual control.

**Returns**
- `phase`: Mean phase calculated over the events cycles.
- `phase_error`: Standard deviation for the mean phase.
"""
function estimate_phase(t, x, events; normalize=true, show_plots=false)

    pr = findpeaks(x, t)
    pks = peakprominences(pr)
    locs = peaklocations(pr)

    if show_plots
        fig, ax = subplots()
        ax.plot(t, x, color="black")
        ax.plot(t[pr.indices], x[pr.indices], "o", color="blue")
        plot_events(events, ax=ax)
    end

    n_events = size(events, 1)
    phase_arr = Float64[]
    for i_event = 2:(n_events-1)
        window_start = events[i_event, 1]
        window_end = events[i_event+1, 1]
        idx = window_start .< locs .< window_end
        window_locs = locs[idx]
        window_pks = pks[idx]
        if isempty(window_pks)
            push!(phase_arr, NaN)
            break
        end
        idx_peak = argmax(window_pks)
        loc = window_locs[idx_peak]
        phase = loc - window_start
        if normalize
            phase /= (window_end - window_start)
        end
        push!(phase_arr, phase)

        if show_plots
            ax.plot(loc, x[pr.indices[idx][idx_peak]], "o", color="red")
        end
    end

    phase = mean(phase_arr)
    phase_error = std(phase_arr)

    return phase, phase_error

end
