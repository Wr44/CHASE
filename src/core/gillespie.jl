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
    new_spacer_chance::Float64,
    θ::Int,
    ε_loop::Float64
)
    t = 0.0 
    tau = 0.0 
    snapshots = Vector{Tuple{Float64, Int, Int, Int, Int}}()
    while t < t_max 
        (rates, lambda) = calculate_rates(bacterias, phages, division_rate, death_rate, infection_rate)

        if (lambda == 0) 
            break 
        end 

        num_bacteria = sum(values(bacterias))
        num_phage = sum(values(phages))

        if min(num_bacteria, num_phage) < θ

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
                            mutation_chance,
                            new_spacer_chance
                        )
                    end
                break
                end 
            end

        else

            tau = (ε_loop * (num_bacteria + num_phage)) / lambda
            
            for (type, args, lambda_i) in rates 
                d = Poisson(lambda_i * tau)
                k = rand(d)

                for _ in 1:k 
                    if type == :division
                        bac_id = args
                        if !haskey(bacterias, bac_id) || bacterias[bac_id] <= 0
                            break
                        end
                        apply_division!(bacterias, bac_id)
                        
                    elseif type == :death
                        bac_id = args
                        if !haskey(bacterias, bac_id) || bacterias[bac_id] <= 0
                            break
                        end
                        apply_death!(bacterias, bac_id)
                        
                    elseif type == :infection_failed
                        bac_id = args[1]
                        phage_id = args[2]
                        if !haskey(bacterias, bac_id) || bacterias[bac_id] <= 0 || 
                        !haskey(phages, phage_id) || phages[phage_id] <= 0
                            break
                        end
                        apply_infection_failed!(bacterias, phages, bac_id, phage_id, new_spacer_chance)
                        
                    elseif type == :infection_succeded
                        bac_id = args[1]
                        phage_id = args[2]
                        if !haskey(bacterias, bac_id) || bacterias[bac_id] <= 0 || 
                        !haskey(phages, phage_id) || phages[phage_id] <= 0
                            break
                        end
                        apply_infection_succeeded!(bacterias, phages, bac_id, phage_id, burst_size, mutation_chance, new_spacer_chance)
                    end
                end
            end 

        end 
    push!(snapshots, (
          t,
          sum(values(bacterias)),
          sum(values(phages)),
          length(bacterias),
          length(phages)
    ))

    t += tau 

    end 

    return snapshots

end