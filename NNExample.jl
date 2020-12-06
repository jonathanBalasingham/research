
include("./Util.jl")
include("./Reaction.jl")
include("./Network.jl")
include("./NetworkModel.jl")
include("./Rates.jl")

using CSV, DifferentialEquations, Flux, CUDA, Optim

reactionsFilepath = "input/reactions_postNN.csv"
speciesFilepath = "input/species_postNN.csv"

reactionsData = CSV.read(reactionsFilepath, DataFrame)
speciesData = CSV.read(speciesFilepath, DataFrame)

ICs = InitialConditions(4.0e-4,2.6e-4,4.6e-4,0)
T=10. 
zeta = 1.
omega = 0.5
F_UV=1.
A_v=2.
E = 0.5
p = Parameters(zeta, omega, T, F_UV, A_v, 0.2, E)

calculateRates!(reactionsData, p)
filterReactionData!(reactionsData, speciesData.name)


network = createNetworkModel(reactionsData, speciesData.name)

u0 = createU0(ICs, network)
tspan = (0.0, 10000000.0)
system = convert(ODESystem, network)


prob = ODEProblem{false}(system, u0, tspan)
@time sol = solve(prob, saveat = 10.)

using Plots

using NeuralPDE
chain = Flux.Chain(Dense(1, 5, σ), Dense(5, 1))
opt = Flux.ADAM(0.1, (0.9, 0.95))
@time sol2 = solve(prob,
                    NNODE(chain, opt), dt=1 / 20f0, verbose=true,
            abstol=1e-10, maxiters=200)
#plot(sol2, yaxis = log10.(sol2.u[0] .+ 1.))
plot(sol, yaxis = log10.(sol.u[1] .+ 1.), size = (1200,800))

