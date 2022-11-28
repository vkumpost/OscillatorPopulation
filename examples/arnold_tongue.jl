## Estimate and plot Arnold tongue
# Load the package
using OscillatorPopulation

# Load model
model = load_model("van-der-pol", "sde")

# Set noise intensity and end simulation time
set_parameter!(model, "σ", 0.005)
set_timespan!(model, 100.0)

# Simulation function
simulation_function = create_simulation_function(
    ["phase_coherence"],
    trajectories=10
)

# Scan a range of input amplitudes and periods
arnold_tongue = scan_arnold(model, simulation_function,
    input_amplitudes=0:0.02:0.4,
    input_periods=0.75:0.02:1.25,
    show_progress=true
)

# Plot the resulting Arnold tongue
plot_arnold(arnold_tongue, "tongue",
    property_name="phase_coherence"
)
