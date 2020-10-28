using DataFrames

function createNetworkModel(reactionsData, species)
    @parameters t r[1:nrow(reactionsData)]
    for s in Symbol.(species) @eval @variables $s(t) end
    reactions = []
    for r in eachrow(reactionsData)
        push!(reactionsData, Reaction(r.rate))
    end
end