using DataFrames
using Catalyst
import SymbolicUtils.FnType

needsSymbol(x) = !isnan(x) && isspecies(x)

# this is heavily relient on the consistency of the reactions csv
function createReaction(row)
    @parameters t
    reactants = Symbol.(filter(needsSymbol,Array(row[1:3])))
    products = Symbol.(filter(needsSymbol, Array(row[4:7])))

    reactant_quantities = repeat([1], length(reactants))
    product_quantities = repeat([1], length(products))

    Catalyst.Reaction(row.rate, 
                      map(x -> (((Sym){(FnType){NTuple{1, Any}, Real}}(x)((ModelingToolkit.value)(t)))), reactants),
                      map(x -> (((Sym){(FnType){NTuple{1, Any}, Real}}(x)((ModelingToolkit.value)(t)))), products),
                      reactant_quantities,
                      product_quantities,
                      only_use_rate = false)
end

function createNetworkModel(reactionsData, species_names)
    @parameters t 
    s = Symbol.(species_names)
    reactions = []
    for r in eachrow(reactionsData)
        push!(reactions, createReaction(r))
    end
    Catalyst.ReactionSystem(reactions,
                            t,
                            map(x -> (((Sym){(FnType){NTuple{1, Any}, Real}}(x)((ModelingToolkit.value)(t)))), s), 
                            [])
end

function createU0(IC::InitialConditions, nw)
    u0 = Float64[]
    for sym in species(nw)
        name = string(sym)
        if name == "C(t)"
            push!(u0, IC.C)
        elseif name == "H(t)"
            push!(u0, IC.H)
        elseif name == "O(t)"
            push!(u0, IC.O)
        elseif name == "N(t)"
            push!(u0, IC.N)
        else
            push!(u0, 0.0)
        end
    end
    u0
end