using Catalyst
using DifferentialEquations


network = @reaction_network begin
    1.69e-35, H2 --> H + H
    6.11271e-61, H + H --> H2
    0.0, H + H2 --> H + H + H2
    0.0, H2 + H2 --> H + H + H2
end

system = convert(ODESystem, network)
u0 = [3.3e-6,0.0]
tspan = (0.0,100.0)
timestamps = collect(range(0,100.0, length=1000))
oprob = ODEProblem(system, u0, tspan)

solve(oprob, TRBDF2(), saveat = timestamps)
