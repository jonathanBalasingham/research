using DataFrames
using Catalyst

needsSymbol(x) = !isnan(x) && isspecies(x)

# this is heavily relient on the consistency of the reactions csv
function createReaction(row)
    @parameters t
    reactants = Symbol.(filter(needsSymbol,Array(row[1:3])))
    products = Symbol.(filter(needsSymbol, Array(row[4:7])))

    reactant_quantities = repeat([1], length(reactants))
    product_quantities = repeat([1], length(products))

    Catalyst.Reaction(row.rate, 
                      map(x -> Catalyst.Variable{Number}(x)(t), reactants),
                      map(x -> Catalyst.Variable{Number}(x)(t), products),
                      reactant_quantities,
                      product_quantities,
                      only_use_rate = false)
end

function createNetworkModel(reactionsData, species_names)
    @parameters t 
    s = Symbol.(species_names)
    #for s in Symbol.(species) @eval @variables $s(t) end
    reactions = []
    for r in eachrow(reactionsData)
        push!(reactions, createReaction(r))
    end
    Catalyst.ReactionSystem(reactions, t, map(x -> Catalyst.Variable{Number}(x)(t), s), [])
end

function createU0(IC::InitialConditions, nw)
    u0 = Float64[]
    for sym in species(nw)
        name = String(sym.name)
        if name == "C"
            push!(u0, IC.C)
        elseif name == "H"
            push!(u0, IC.H)
        elseif name == "O"
            push!(u0, IC.O)
        else
            push!(u0, 0.0)
        end
    end
    u0
end