## Simulate a population of uncoupled stochastic oscillators
# Load the package
using OscillatorPopulation

# Load model and change some of its parameters
model = load_model("van-der-pol", "sde")
set_parameter!(model, ["I", "σ"], [0.5, 0.08])

# Generate input signal (events) and add it to the model
events = create_events(
    [(:DD, 5), (:LD, 5, 0.5, 0.5), (:LL, 5)]
)
set_input!(model, events)

# Set maximum integration time to the end of the input signal
set_timespan!(model, events[end])

# Simulate population of 1000 oscillators
solution = simulate_population(model, 1000, seed=3)

# Plot solution
plot_solution(solution)
