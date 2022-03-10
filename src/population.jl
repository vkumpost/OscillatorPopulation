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
`simulate_population(model::Model, trajectories=1; save_trajectories=true,
    show_progress=false)`

Simulate a population.

**Arguments**
- `model`: Model.

**Optional Arguments**
- `trajectories`: Number of trajectories to simulate.

**Keyword Arguments**
- `save_trajectories`: If `false`, only mean is saved.
- `show_progress`: If `true`, show a progress bar in the terminal.
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
