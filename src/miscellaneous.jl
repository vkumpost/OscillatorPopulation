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