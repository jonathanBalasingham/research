include("./Util.jl")
include("./Reaction.jl")
include("./Network.jl")
#include("./NetworkModel.jl")
include("./Rates.jl")

using CSV
using DifferentialEquations

function formulate_problem_plain_ODE(reactionsFilepath::String, speciesFilepath::String)
    println(reactionsFilepath)
    println(speciesFilepath)
    reactionsData = CSV.read(reactionsFilepath)
    speciesData = CSV.read(speciesFilepath)

    # These will be an input eventually
    ICs = InitialConditions(1.0,2.6e-4,4.6e-4)
    T=10 
    zeta = 1.3e-17
    F_UV=1
    A_v=10
    p = Parameters(zeta, 1.0, 1.0, T, F_UV, A_v, 1.0, 0.5)

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
    p = Parameters(zeta, 1.0, 1.0, T, F_UV, A_v, 1.0, 0.5)

    calculateRates!(reactionsData, p)
    filterReactionData!(reactionsData, speciesData["name"])
    network = createNetworkModel(reactionsData, speciesData["name"])
    u0 = createU0(ICs, network)
    tspan = (0.0, 10.0)
    system = convert(ODESystem, network)
    ODEProblem(system, u0, tspan)
end

#julia_main(ARGS[1], ARGS[2])
