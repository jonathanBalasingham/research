

struct Species
    name::String
    mass::Int
    natoms::Int
end

function createSpeciesList(speciesData)::Array{Species}
    speciesList = []
    for i in eachrow(speciesData)
        push!(speciesList, Species(i.name, i.mass, i.natoms))
    end
    speciesList
end