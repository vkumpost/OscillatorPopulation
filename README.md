# OscillatorPopulation

Software for my thesis.

## Installation

The package can be installed by running
```julia
import Pkg
Pkg.add(url="https://github.com/vkumpost/OscillatorPopulation")
```

To make sure everything is ready to go we can run package tests
```julia
Pkg.test("OscillatorPopulation")
```

If the compilation fails due to `PyCall` package not being able to find a Python
installation, a simple fix might be to run
```julia
ENV["PYTHON"]=""
Pkg.build("PyCall")
```
