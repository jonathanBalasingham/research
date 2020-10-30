using DataFrames
using Catalyst

quote
    begin
        var"#570#t" = (Catalyst.Variable{ModelingToolkit.Parameter{Number}}(:t))()
        (var"#570#t",)
    end
    begin
        var"#571#X" = (Catalyst.Variable{Number}(:X))(var"#570#t")
        var"#572#Y" = (Catalyst.Variable{Number}(:Y))(var"#570#t")
        var"#573#XY" = (Catalyst.Variable{Number}(:XY))(var"#570#t")
        (var"#571#X", var"#572#Y", var"#573#XY")
    end
    Catalyst.ReactionSystem([Catalyst.Reaction(1.0, [var"#571#X", var"#572#Y"], [var"#573#XY"], [1, 1], [1], only_use_rate = false)], var"#570#t", [var"#571#X", var"#572#Y", var"#573#XY"], [])
end

makeSymbol(x) = !isnan(x) && isspecies(x)

function createReaction(row)
    reactants = Symbol.(filter(makeSymbol,Array(row[1:3])))
    reactant_quantities = repeat([1], length(reactants))
    products = Symbol.(filter(makeSymbol, Array(row[4:7])))
    product_quantities = repeat([1], length(products))
    Catalyst.Reaction(row.rate, 
                      Catalyst.Variable{Number}.(reactants)(t),
                      Catalyst.Variable{Number}.(products)(t),
                      reactant_quantities,
                      product_quantities,
                      only_use_rate = false)
end

function createNetworkModel(reactionsData, species)
    @parameters t r[1:nrow(reactionsData)]
    for s in Symbol.(species) @eval @variables $s(t) end
    reactions = []
    for r in eachrow(reactionsData)
        push!(reactions, createReaction(r))
    end
end