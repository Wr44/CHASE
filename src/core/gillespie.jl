using Distributions

function gillespie(
    bacterias::Bacterias,
    phages::Phages,
    t_max::Int,
    division_rate::Float64,
    death_rate::Float64,
    burst_size::Int,
    infection_rate::Float64,
    mutation_chance::Float64,
    new_spacer_chance::Float64
)
    t = 0.0 
    snapshots = Vector{Tuple{Float64, Int, Int}}() 
    while t < t_max 
        (rates, lambda) = calculate_rates(bacterias, phages, division_rate, death_rate, infection_rate)
        
        if (lambda == 0) 
            break 
        end 

        probas = cumsum([x[3]/lambda for x in rates])

        d = Exponential(1.0/lambda)
        tau = rand(d)

        r = rand()

        for i in eachindex(probas) 
            if r < probas[i] 
                if (rates[i][1] == :division)
                    apply_division!(
                        bacterias,
                        rates[i][2]
                    )
                elseif (rates[i][1] == :death)
                    apply_death!(
                        bacterias,
                        rates[i][2]                        
                    )
                elseif (rates[i][1] == :infection_failed)
                    apply_infection_failed!(
                            bacterias, 
                            phages, 
                            rates[i][2][1], 
                            rates[i][2][2], 
                            new_spacer_chance
                    )
                elseif (rates[i][1] == :infection_succeded)
                    apply_infection_succeeded!(
                            bacterias, 
                            phages, 
                            rates[i][2][1], 
                            rates[i][2][2], 
                            burst_size, 
                            mutation_chance
                    )
                end
                break
            end 
        end

        push!(snapshots, (t, sum(values(bacterias)), sum(values(phages))))
        t += tau 
    end

    return snapshots

end