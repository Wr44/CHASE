module CHASE

export Bacterias, Phages, recognize
export apply_division!, apply_death!, apply_phage_decay!
export apply_infection_succeeded!, apply_infection_failed!
export init_rates, update_rates, gillespie

include("core/state.jl")
include("core/rates.jl")
include("core/events.jl")
include("core/gillespie.jl")

end
