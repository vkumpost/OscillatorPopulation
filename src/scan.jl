"""
`_find_all_combinations`

Find all possible combinations among the vectors of values. 

**Arguments**
- `vectors`: Vector of vectors representing different value sets.

**Returns**
- `combinations`: Matrix containing all combinations of values from the
    individual vectors. Each row represents one combination.

**Example**
```
julia> vectors = [[1, 2], [3, 4, 5], [6]];
julia> combinations = _find_all_combinations(vectors);
julia> print(combinations);
[1 3 6; 1 4 6; 1 5 6; 2 3 6; 2 4 6; 2 5 6]
```
"""
function _find_all_combinations(vectors)

    # Number of input vectors
    n = length(vectors)  

    # The last column of the combinations matrix is the last input vector
    combinations = hcat(vectors[n])

    # Iterate the remaining vectors
    for i = (n-1):-1:1

        # Get the current vector and its length
        x = vectors[i]
        nx = length(x)

        # Get the current size of the combinations matrix
        nc = size(combinations, 1)

        # Repeat the current combinations for each element in x 
        combinations = repeat(combinations, nx)

        # Add an element from x to each repetition of combinations
        x = vcat([fill(xi, nc) for xi in x]...)
        combinations = hcat(x, combinations)

    end

    return combinations

end


"""
`scan`

Scan parameters of a model.

**Arguments**
- `model`: Model.
- `parameters`: A vector of pairs
- `simulation_function`: A function that takes in `PopulationSolution`, and returns
    some descriptive parameters. If the function is called without any arguments,
    it returns the names of the parameters as a vector of strings. For example:
    ```
    julia> simulation_function(solution)
    [1.1, 2.2, 3.3]
    julia> simulation_function()
    ["Period", "Amplitude", "Phase"]
    ```

**Keyword Arguments**
- `show_progress`: If `true`, show a progress bar in the terminal.

**Returns**
- `df`: Dataframe with results. First columns correspond to the input parameters
    followed by columns obtained by `summary_function` for the corresponding
    parameters.
"""
function scan(model::Model, parameters::Vector, simulation_function::Function;
    show_progress::Bool=false)

    # Get values from the parameter dictionary
    parameter_names = [x[1] for x in parameters]
    parameter_values = [x[2] for x in parameters]

    # Find all possible parameter combinations
    parameter_value_combinations = _find_all_combinations(parameter_values)
    n = size(parameter_value_combinations, 1)
    
    # Initialize summary matrix
    summary_names = simulation_function()
    scan_summary = fill(NaN, n, length(summary_names))

    # Initialize progress bar
    if show_progress
        progressmeter = ProgressMeter.Progress(n; barlen=20)
    end

    # Iterate all parameter value combinations
    lk = ReentrantLock()
    Threads.@threads for i = 1:n

        # Protect the original model from overwriting
        model2 = deepcopy(model)

        # Set model parameters
        set_parameter!(model2, parameter_names, parameter_value_combinations[i, :])

        # Evaluate the sumary function
        simulation_summary = simulation_function(model2)

        # Save the model summary to the overall matrix
        scan_summary[i, :] = simulation_summary

        # Update progress bar
        if show_progress
            lock(lk) do
                ProgressMeter.next!(progressmeter)
            end
        end

    end

    # Build the output dataframe
    matrix = hcat(parameter_value_combinations, scan_summary)
    names = vcat(parameter_names, summary_names)
    df = DataFrame(matrix, names)

    return df

end
