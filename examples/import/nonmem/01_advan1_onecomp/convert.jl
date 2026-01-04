# NONMEM ADVAN1 Import - Julia Example
# Run: julia --project=core/OpenPKPDCore convert.jl

using OpenPKPDCore
using JSON

println("NONMEM ADVAN1 Import")
println("="^50)

# Import NONMEM control file
model = import_nonmem("run001.ctl")

println("\nImported Model:")
println("  Type:       $(model.model_type)")
println("  Parameters: $(model.parameters)")

# Display OpenPKPD specification
println("\nOpenPKPD Specification:")
println("-"^50)
spec = model.spec

println("  Model:   $(spec.model)")
println("  Params:  $(spec.params)")

if !isnothing(spec.iiv)
    println("  IIV:")
    println("    Kind:   $(spec.iiv.kind)")
    println("    Omegas: $(spec.iiv.omegas)")
end

if !isnothing(spec.error)
    println("  Error:   $(spec.error)")
end

# Compare with expected
println("\n" * "="^50)
println("VALIDATION")
println("="^50)

expected = JSON.parsefile("expected.json")

# Validate parameters
params_match = spec.params == expected["params"]
println("Parameters match: $params_match")

# Validate IIV
iiv_match = spec.iiv.omegas == expected["iiv"]["omegas"]
println("IIV match: $iiv_match")

# Simulate to verify
println("\n" * "="^50)
println("SIMULATION TEST")
println("="^50)

grid = SimulationGrid(t0=0.0, t1=24.0, saveat=0:1.0:24)
result = simulate(spec, grid)

println("\nSimulation successful: $(length(result.times)) time points")
println("Cmax: $(round(maximum(result.observations["conc"]), digits=3)) mg/L")
