## Find parameter values that give 24-hour oscillations
# Load the package
using OscillatorPopulation

# Load model
model = load_model("kim-forger-full", "ode")
set_output!(model, 2)

# Specify parameters and their search range
parameter_names = ["A", "dM=dP=dR", "I"]
search_range = [(0.0, 100.0), (0.0, 1.0), (0.0, 0.1)]

# Target period estimated from data
target_period = 24.0

# Construct cost function
cost_function = create_entrainable_oscillator_objective(model, target_period,
    parameter_names=parameter_names
)

# Perform optimization
best_candidate, final_population = optimize(cost_function,
    search_range=search_range,
    max_steps=5_000
)
