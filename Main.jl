include("./Util.jl")
include("./Reaction.jl")
include("./Network.jl")
include("./NetworkModel.jl")
include("./Rates.jl")

using CSV
using DifferentialEquations
using Plots


function formulate_problem_plain_ODE(reactionsFilepath::String, speciesFilepath::String)
    reactionsData = CSV.read(reactionsFilepath)
    speciesData = CSV.read(speciesFilepath)

    # These will be an input eventually
    ICs = InitialConditions(1.0,2.6e-4,4.6e-4)
    T=10 
    zeta = 1.3e-17
    F_UV=1
    A_v=10
    p = Parameters(zeta, 0.4, 0.3, T, F_UV, A_v, 0.2, 0.5)

    calculateRates!(reactionsData, p)
    filterReactionData!(reactionsData, speciesData["name"])
    prob = createNetwork(ICs, speciesData["name"], reactionsData)
    #sol = solve(prob)
end

function formulate_problem_catalyst(rfp::String, sfp::String)
    reactionsData = CSV.read(rfp)
    speciesData = CSV.read(sfp)

    ICs = InitialConditions(1.0,2.6e-4,4.6e-4)
    T=10 
    zeta = 1.3e-17
    F_UV=1
    A_v=10
    p = Parameters(zeta, 0.4, 0.3, T, F_UV, A_v, 0.2, 0.5)

    calculateRates!(reactionsData, p)
    filterReactionData!(reactionsData, speciesData["name"])
    network = createNetworkModel(reactionsData, speciesData["name"])
    u0 = createU0(ICs, network)
    tspan = (0.0, 10000000.0)
    system = convert(ODESystem, network)
    ODEProblem(system, u0, tspan)
end

#julia_main(ARGS[1], ARGS[2])

# speed test
rfp = "input/reactions.csv"
sfp = "input/species.csv"

sa = collect(range(0.0, 10000000.0, step = 0.1))

@time sol1 = solve(formulate_problem_catalyst(rfp, sfp), saveat=10.)
#@time sol2 = solve(formulate_problem_plain_ODE(rfp,sfp), saveat=sa)

plot(sol1)