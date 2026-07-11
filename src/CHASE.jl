module CHASE

export Bacterias, Phages, recognize
export apply_division!, apply_death!
export apply_infection_succeeded!, apply_infection_failed!
export calculate_rates, gillespie

include("core/state.jl")
include("core/rates.jl")
include("core/events.jl")
include("core/gillespie.jl")

end
