include("./Species.jl")
include("./Reaction.jl")

using CSV


function julia_main(reactionsFilepath::AbstractString, speciesFilepath::AbstractString)
    try
        reactionsData = CSV.read(reactionsFilepath)
        speciesData = CSV.read(speciesFilepath)
    catch e
        println("Error occured trying to read input files")
        println(e)
    end

    # These will be an input eventually
    ICs = InitialConditions(1.0,2.6e-4,4.6e-4)
    T=10 
    zeta = 1.3e-17
    F_UV=1
    A_v=10
    p = Parameters(zeta,1.0,1.0,T,F_UV,A_v, 1.0, 0.5)

    calculateRates!(reactionsData, p)
    filterReactionData!(reactionsData, speciesData["name"])
    res = createNetwork(ICs, speciesData["name"], reactionsData)

    speciesList = createSpeciesList(speciesData)
    reactionList = createReactionList(reactionsData, speciesData["name"])
end
