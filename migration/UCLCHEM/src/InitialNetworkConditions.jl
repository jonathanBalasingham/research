
struct InitialConditions
    H::Float64
    C::Float64
    O::Float64
    N::Float64
end


struct InitialNetworkConditions
    initialConcentrations::Dict{String, Float64}
    settingsFilepath::String
end