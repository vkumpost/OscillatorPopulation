"""
`smooth`

Smooth a vector using a moving average filter. Endpoints are handled by
collapsing the length of the filter as showed below.
```
yy[1] = y[1]
yy[2] = (y[1] + y[2] + y[3]) / 3
yy[3] = (y[1] + y[2] + y[3] + y[4] + y[5]) / 5
yy[4] = (y[2] + y[3] + y[4] + y[4] + y[6]) / 5
```

**Arguments**
- `y`: Array of values.

**Keyword Arguments**
- `span`: Filter length (default: `5`).

**Returns**
- `yy`: Smoothed array of values.
"""
function smooth(y; span=5)

    half_span = floor(Int, span/2)
    n = length(y)
    yy = fill(NaN, n)

    for i = 1:n

        i_span = min(i - 1, n - i, half_span)
        i_left = i - i_span
        i_right = i + i_span

        yy[i] = mean(y[i_left:i_right])

    end

    return yy

end


"""
`generate_random_values`

Generate a vector of random numbers drawn from the truncated normal distribution.

**Arguments**
- `μ`: Mean.
- `σ`: Standard deviation.
- `n`: Number of values.

**Keyword Arguments**
- `lower`: Minimum value (default: `nothing`).
- `upper`: Maximum value (default: `nothing`).
- `seed`: Random number generator seed (default: `nothing`).

**Returns**
- `values`: An array of generated values.
"""
function generate_random_values(μ, σ, n; lower=nothing, upper=nothing, seed=nothing)

    random_number_generator = MersenneTwister(seed)
    distribution = truncated(Normal(μ, σ), lower, upper)
    values = rand(random_number_generator, distribution, n)

    return values
    
end


"""
`benchmark`

Estimate evaluation time for a function.

**Arguments**
- `fun`: Function.

**Keyword Arguments**
- `n_repeats`: Number of repeats (default: 100).

**Returns**
- `mean_time`: Mean evaluation time in milliseconds.
- `std_time`: Standard deviation in milliseconds.
"""
function benchmark(fun; n_repeats=100)
    
    # Dummy run to make sure the function is compiled
    _ = @elapsed x = fun()

    time_array = fill(NaN, n_repeats)
    for i = 1:n_repeats
        time_array[i] = @elapsed x = fun()
    end

    mean_time = mean(time_array .* 1000)
    std_time = std(time_array .* 1000)

    return mean_time, std_time

end


"""
`find_closest`

Find the index of the closest value in an array.

**Arguments**
- `array`: Array of values.
- `value`: Target value.

**Returns**
- `index`: Index of the value in `array` that is the closest to `value`.
"""
function find_closest(array, value)

    difference = abs.(array .- value)
    index = argmin(difference)
    return index
    
end


"""
`binary_search`

Use binary search to find `x` for which `fun` evaluates to `target_value`.

**Arguments**
- `fun`: A monotonically increasing function which takes one number as a parameter.
- `target_value`: Targer value to which should `fun` evaluate.
- `search_range`: A two-element array specifying the search range.

**Keyword Arguments**
- `tolerance`: Maximal tolerance from the target value.
- `max_steps`: Maximal number of steps.
- `verbose`: If `true`, print the estimated value and error at each step.

**Returns**
- `x`: A number for which `fun(x) ≈ target_value`.
"""
function binary_search(fun, target_value, search_range; tolerance=√eps(),
    max_steps=1_000_000, verbose=false)
    
    # Get search boundaries
    x_left = search_range[1]
    x_right = search_range[2]

    # Check that the target value lies in the search interval
    if !(fun(x_left) < target_value < fun(x_right))
        msg = "Target value is not in the search interval!"
        err = OscillatorPopulationError(msg)
        throw(err)
    end

    # Initialized loop variables
    estimated_value = fun(x_left)
    x_new = x_left
    step_counter = 0

    while (abs(estimated_value - target_value) > tolerance) && (step_counter < max_steps)
        
        # Guess new value in the middle of the search interval
        x_new = (x_left + x_right) / 2
        estimated_value = fun(x_new)

        # Adjust the search interval
        if estimated_value < target_value
            x_left = x_new
        else
            x_right = x_new
        end

        # Increase step counter
        step_counter += 1

        # Print step number, estimated value and error
        if verbose
            println("Step $(step_counter), Estimate $(round(x_new, digits=10)), Error $(round(abs(estimated_value - target_value), digits=10))")
        end

    end

    return x_new

end
