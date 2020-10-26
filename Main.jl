include("./Util.jl")
include("./Species.jl")
include("./Reaction.jl")
include("./Network.jl")
include("./Rates.jl")
include("./Chemistry.jl")

using CSV
using DifferentialEquations

function julia_main(reactionsFilepath::String, speciesFilepath::String)
    println(reactionsFilepath)
    println(speciesFilepath)
    try
        reactionsData = CSV.read(reactionsFilepath)
        speciesData = CSV.read(speciesFilepath)
    catch e
        println("Error occured trying to read input files")
        println(e)
    end
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
    sol = solve(prob)
end

julia_main(ARGS[1], ARGS[2])
