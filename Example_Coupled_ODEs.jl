using DifferentialEquations, ModelingToolkit

@parameters t
@variables y(t) 
@variables x(t)
@derivatives D'~t
eqs = [
    D(y) ~ -y + 3*x,
    D(x) ~ 4*x-2*y
]

coupled = ODESystem(eqs, name=:coupled)

u0 = [x => 1.0,
      y => 2.0]
tspan = (0.0,1.0)

prob = ODEProblem(coupled, u0, tspan)
solve(prob)
