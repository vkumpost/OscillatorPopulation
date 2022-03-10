"""
`PopulationSolution`

Solution for a population of models.

**Fields**
- `time`: Time vector (time steps).
- `mean`: Population mean over individual trajectories (time steps x variables).
- `trajectories`: Individual trajectories (time steps x variables x trajectories).
- `success`: True, if the integration was successful.
"""
struct PopulationSolution

    time::Vector
    mean::Matrix
    trajectories::Array
    events::Matrix
    success::Bool

end


"""
`simulate_population`

Simulate a population.

**Arguments**
- `model`: Model.

**Optional Arguments**
- `trajectories`: Number of trajectories to simulate.

**Keyword Arguments**
- `save_trajectories`: If `false`, only mean is saved.
- `show_progress`: If `true`, show a progress bar in the terminal.

**Returns**
- `solution`: PopulationSolution.
"""
function simulate_population(model::Model, trajectories=1;
    save_trajectories=true, show_progress=false)

    # Prepare variables for PopulationSolution fields
    t = Vector{Float64}(undef, 0)  # time
    m = Matrix{Float64}(undef, 0, 0)  # mean
    U = Array{Float64, 3}(undef, 0, 0, 0)  # trajectories
    events = model.input[1]  # events
    
    # Initialize the porgress meter
    if show_progress
        progress_meter = ProgressMeter.Progress(trajectories; barlen=20)
    end

    # Simulate trajectories
    success_arr = fill(false, trajectories)  # success for individual trajectories
    lk = ReentrantLock()
    Threads.@threads for i = 1:trajectories

        # Make a copy of the model
        model2 = deepcopy(model)

        # Solve the problem
        sol_x, sol_t, sol = simulate_model(model2)

        # Check if the simulation was successful
        if model2.problem isa Union{ODEProblem, SDEProblem}
            success_arr[i] = sol.retcode == :Success
        elseif model2.problem isa JumpProblem
            success_arr[i] = sol.retcode == :Default
        else
            error("Unknown model type!")
        end

        # Save the solution if the current simulation was successful
        if success_arr[i]

            lock(lk) do

                # Save time
                if isempty(t)
                    t = sol_t
                end

                # Add the trajectory contribution to the mean
                if isempty(m)
                    # Initialize the mean vector
                    m = (sol_x ./ trajectories)
                else
                    # Add solution to the mean vector
                    m .+= (sol_x ./ trajectories)
                end

                # Save the individual trajectory
                if save_trajectories
                    if isempty(U)
                        # Initialize the matrix for individual trajectories
                        nsamples = size(sol_x, 1)
                        nvariables = size(sol_x, 2)
                        U = fill(NaN, nsamples, nvariables, trajectories)
                    end
                    # Add the current trajectory to the trajectory list
                    U[:, :, i] = sol_x
                end

            end

        else
            # Stop the simulation if the current simulation was unsuccessful
            break
        end

        if show_progress
            lock(lk) do
                ProgressMeter.next!(progress_meter)
            end
        end

    end

    # Build the output struct
    success = all(success_arr)
    solution = PopulationSolution(t, m, U, events, success)

    return solution

end


"""
`plot_solution`

Plot a solution.

**Arguments**
- `solution`: PopulationSolution.

**Optional Arguments**
- `variable`: Index of the variable to be plotted.

**Keyword Arguments**
- `n`: Number of trajectories to plot in the background.
- `ax`: PyPlot axes.

**Returns**
- `ax`: PyPlot axes.
"""
function plot_solution(solution::PopulationSolution, variable=1; n=3, ax=gca())

    t = solution.time
    m = solution.mean
    U = solution.trajectories
    n_trajectories = size(U, 3)

    if !isempty(U)
        for i = 1:min(n, n_trajectories)
            x = U[:, variable, i]
            ax.plot(t, x, color="gray", alpha=0.5)
        end
    end

    x = m[:, variable]
    ax.plot(t, x, color="black")

    plot_events(solution.events, ax=ax)

    ax.set_xlabel("Time")
    ax.set_ylabel("Variable $(variable)")

    return ax

end


"""
`select_time`

Select a specific time from a solution.

**Arguments**
- `solution`: `PopulationSolution`.

**Keyword Arguments**
- `remove_offset`: If `true`, the first time point is set to 0.
- `offset`: Set offset (the first time point).
- `min_time`: Minimal time (inclusive).
- `max_time`: Minimal time (exclusive).

**Returns**
- `PopulationSolution`
"""
function select_time(solution::PopulationSolution; remove_offset=true, kwargs...)

    solution = deepcopy(solution)

    t = solution.time
    m = solution.mean
    U = solution.trajectories
    events = solution.events
    success = solution.success

    if (:min_time in keys(kwargs))
        
        # Remove samples before `min_time`
        min_time = kwargs[:min_time]
        indices = t .>= min_time
        t = t[indices]
        m = m[indices, :]
        if !isempty(U)
            U = U[indices, :, :]
        end

        # Remove events before `min_time`
        new_events = Matrix{Float64}(undef, 0, 2)
        for i = 1:size(events, 1)
            if events[i, 2] >= min_time
                if events[i, 1] >= min_time
                    new_events = vcat(new_events, events[i, :]')
                elseif min_time != events[i, 2]
                    new_events = vcat(new_events, [min_time events[i, 2]])
                end
            end 
        end
        events = new_events

    end

    if (:max_time in keys(kwargs))

        # Remove samples after `max_time`
        max_time = kwargs[:max_time]
        indices = t .< max_time
        t = t[indices]
        m = m[indices, :]
        if !isempty(U)
            U = U[indices, :, :]
        end

        # Remove events after `max_time`
        newevents = Matrix{Float64}(undef, 0, 2)
        for i = 1:size(events, 1)
            if events[i, 1] < max_time
                if events[i, 2] < max_time
                    newevents = vcat(newevents, events[i, :]')
                elseif events[i, 1] != max_time
                    newevents = vcat(newevents, [events[i, 1] max_time])
                end
            end 
        end
        events = newevents
    end

    if remove_offset
        # Set first timepoint to zero
        events .-= t[1]
        t .-= t[1]
    end

    if :offset in keys(kwargs)
        # Set first timepoint to offset
        offset = kwargs[:offset]
        events .+= offset .- t[1]
        t .+= offset .- t[1]
    end

    return PopulationSolution(t, m, U, events, success)

end
