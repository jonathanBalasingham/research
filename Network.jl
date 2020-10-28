using ModelingToolkit
# what is D??
function prod(ind, species, reactions, Y)
    name = species[ind]
    # find reactions that have name in prod 1,2,3 or 4
    rows = filter(row -> name in row[4:7], reactions)
    p = 0
    for r in eachrow(rows)
        re1,re2,re3 = r[1:3]
        term = 1
        for (i,s) in enumerate(species)
            if s == re1
                term *= Y[i]
            elseif s == re2
                term *= Y[i]
            elseif s == re3
                term *= Y[i]
            end
        end
        p += term*r.rate
    end
    p
end

function loss(ind, species, reactions, Y)
    name = species[ind]
    # find reactions that have name in prod 1,2,3 or 4
    rows = filter(row -> name in row[1:3], reactions)
    # filter the row
    p = 0
    for r in eachrow(rows)
        pr1,pr2,pr3,pr4 = r[4:7]
        term = -1
        for (i,s) in enumerate(species)
            if s == pr1
                term *= Y[i]
            elseif s == pr2
                term *= Y[i]
            elseif s == pr3
                term *= Y[i]
            elseif s == pr4
                term *= Y[i]
            end
        end
        p += term*r.rate
    end
    p
end

function createNetwork(IC::InitialConditions, species::Array{String}, reactionsData)
    ModelingToolkit.@parameters t
    @variables y[1:length(species)](t)
    @derivatives D'~t

    ydot = D.(y)
    eqs = []
    u0 = []
    for (i,yt) in enumerate(y)
        name = species[i]
        if name == "C"
            push!(u0, yt => IC.C)
        elseif name == "H"
            push!(u0, yt => IC.H)
        elseif name == "O"
            push!(u0, yt => IC.O)
        else
            push!(u0, yt => 0.0)
        end
        push!(eqs, ydot[i] ~ prod(i,species,reactionsData, y) + y[i]*loss(i,species,reactionsData,y))
    end
    eqs = eqs .|> simplify
    println(eqs)
    tspan = (0.0, 1.0)
    network = ODESystem(eqs, name=:uclchem)
    ODEProblem(network, u0, tspan)
end