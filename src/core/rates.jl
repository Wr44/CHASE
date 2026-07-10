function calculate_rates(
    bacterias::Bacterias,
    phages::Phages,
    division_rate::Float64,
    death_rate::Float64,
    infection_rate::Float64,
)::Tuple{Vector{Tuple{Symbol, Any, Float64}}, Float64}

    v = Vector{Tuple{Symbol, Any, Float64}}()
    
    for (spacers, bacteria_count) in bacterias
        push!(v, (:division, spacers, division_rate * bacteria_count))
        push!(v, (:death, spacers, death_rate * bacteria_count))

        for (phage_id, phage_count) in phages
            if recognize(spacers, phage_id)
                push!(v, (:infection_failed, (phage_id, spacers) , infection_rate*bacteria_count*phage_count))
            else
                push!(v, (:infection_succeded, (phage_id, spacers) , infection_rate*bacteria_count*phage_count))
            end    
        end
    end

    lambda = sum(t -> t[3], v)
    return (v, lambda)
end