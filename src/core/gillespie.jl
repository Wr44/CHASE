using Distributions

function gillespie(
    bacterias::Bacterias,
    phages::Phages,
    t_max::Int,
    division_rate::Float64,
    death_rate::Float64,
    phage_decay::Float64,
    K::Int,
    burst_size::Int,
    infection_rate::Float64,
    mutation_chance::Float64,
    new_spacer_chance::Float64,
    θ::Int,
    ε_loop::Float64
)
    t = 0.0
    tau = 0.1
    snapshots = Vector{Tuple{Float64, Int, Int, Int, Int}}()

    cache = init_rates(bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate)

    while t < t_max

        if cache.lambda == 0
            break
        end

        num_bacteria = sum(values(bacterias))
        num_phage = sum(values(phages))

        if min(num_bacteria, num_phage) < θ

            all_rates = vcat(
                [(:division, k, v) for (k, v) in cache.division],
                [(:death, k, v) for (k, v) in cache.death],
                [(:phage_decay, k, v) for (k, v) in cache.phage_decay],
                [(:infection_failed, k, v) for (k, v) in cache.infection_failed],
                [(:infection_succeeded, k, v) for (k, v) in cache.infection_succeeded]
            )

            probas = cumsum([x[3] / cache.lambda for x in all_rates])
            tau = rand(Exponential(1.0 / cache.lambda))
            r = rand()

            for i in eachindex(probas)
                if r < probas[i]
                    event_type = all_rates[i][1]
                    event_data = all_rates[i][2]

                    if event_type == :division
                        apply_division!(bacterias, event_data)
                        update_rates!(cache, bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate, :division, event_data)

                    elseif event_type == :death
                        apply_death!(bacterias, event_data)
                        update_rates!(cache, bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate, :death, event_data)

                    elseif event_type == :phage_decay
                        apply_phage_decay!(phages, event_data)
                        update_rates!(cache, bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate, :phage_decay, event_data)

                    elseif event_type == :infection_failed
                        spacers, phage_id = event_data
                        apply_infection_failed!(bacterias, phages, spacers, phage_id, new_spacer_chance)
                        update_rates!(cache, bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate, :phage_decay, phage_id)

                    elseif event_type == :infection_succeeded
                        spacers, phage_id = event_data
                        old_bac_keys = Set(keys(bacterias))
                        old_phage_keys = Set(keys(phages))

                        apply_infection_succeeded!(bacterias, phages, spacers, phage_id, burst_size, mutation_chance, new_spacer_chance)

                        update_rates!(cache, bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate, :death, spacers)
                        update_rates!(cache, bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate, :phage_decay, phage_id)

                        for new_spacers in setdiff(keys(bacterias), old_bac_keys)
                            update_rates!(cache, bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate, :new_clone, new_spacers)
                        end

                        for new_phage_id in setdiff(keys(phages), old_phage_keys)
                            update_rates!(cache, bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate, :new_phage, new_phage_id)
                        end
                    end

                    break
                end
            end

        else

            tau_calc = (ε_loop * (num_bacteria + num_phage)) / cache.lambda
            tau_max = ε_loop / division_rate
            tau = min(tau_calc, tau_max)

            all_rates = vcat(
                [(:division, k, v) for (k, v) in cache.division],
                [(:death, k, v) for (k, v) in cache.death],
                [(:phage_decay, k, v) for (k, v) in cache.phage_decay],
                [(:infection_failed, k, v) for (k, v) in cache.infection_failed],
                [(:infection_succeeded, k, v) for (k, v) in cache.infection_succeeded]
            )

            rates_with_k = map(r -> (r[1], r[2], rand(Poisson(r[3] * tau))), all_rates)

            for (type, args, k) in rates_with_k
                for _ in 1:k
                    if type == :death
                        if !haskey(bacterias, args) || bacterias[args] <= 0
                            break
                        end
                        apply_death!(bacterias, args)

                    elseif type == :division
                        if !haskey(bacterias, args) || bacterias[args] <= 0
                            break
                        end
                        apply_division!(bacterias, args)

                    elseif type == :phage_decay
                        if !haskey(phages, args) || phages[args] <= 0
                            break
                        end
                        apply_phage_decay!(phages, args)

                    elseif type == :infection_failed
                        bac_id, phage_id = args
                        if !haskey(bacterias, bac_id) || bacterias[bac_id] <= 0 ||
                           !haskey(phages, phage_id) || phages[phage_id] <= 0
                            break
                        end
                        apply_infection_failed!(bacterias, phages, bac_id, phage_id, new_spacer_chance)

                    elseif type == :infection_succeeded
                        bac_id, phage_id = args
                        if !haskey(bacterias, bac_id) || bacterias[bac_id] <= 0 ||
                           !haskey(phages, phage_id) || phages[phage_id] <= 0
                            break
                        end
                        apply_infection_succeeded!(bacterias, phages, bac_id, phage_id, burst_size, mutation_chance, new_spacer_chance)
                    end
                end
            end

            cache = init_rates(bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate)
        end

        t += tau

        push!(snapshots, (
            t,
            sum(values(bacterias)),
            sum(values(phages)),
            length(bacterias),
            length(phages)
        ))
    end

    return snapshots
end