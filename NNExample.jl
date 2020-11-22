
include("./Util.jl")
include("./Reaction.jl")
include("./Network.jl")
include("./NetworkModel.jl")
include("./Rates.jl")

reactionsFilepath = "input/reactions.csv"
speciesFilepath = "input/species.csv"

reactionsData = CSV.read(reactionsFilepath)
speciesData = CSV.read(speciesFilepath)

ICs = InitialConditions(4.0e-4,2.6e-4,4.6e-4,0)
T=1. 
zeta = 1.3e-17
F_UV=1.
A_v=1.
p = Parameters(zeta, 0.4, 0.3, T, F_UV, A_v, 0.2, 0.5)

calculateRates!(reactionsData, p)
filterReactionData!(reactionsData, speciesData["name"])


network = createNetworkModel(reactionsData, speciesData["name"])

u0 = createU0(ICs, network)
tspan = (0.0, 10000.0)
system = convert(ODESystem, network)


prob = ODEProblem(system, u0, tspan)
sol = solve(prob, saveat = 10.)

using Plots
plot(sol)

using NeuralPDE
chain = Flux.Chain(Dense(1, 5, σ), Dense(5, 1))
opt = Flux.ADAM(0.1, (0.9, 0.95))
@time sol2 = solve(prob,
                    NeuralPDE.NNRODE(chain, opt), dt=1 / 20f0, verbose=true,
            abstol=1e-10, maxiters=200)
plot(sol2)

