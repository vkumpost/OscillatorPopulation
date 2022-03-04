"""
`OscillatorPopulationError`

A generic exception for this package.
"""
struct OscillatorPopulationError <: Exception
    msg::String
end

function Base.showerror(io::IO, err::OscillatorPopulationError)
    print(io, "OscillatorPopulationError: ")
    print(io, err.msg)
end
