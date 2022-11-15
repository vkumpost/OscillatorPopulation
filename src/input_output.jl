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
