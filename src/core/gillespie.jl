using Distributions

function gillespie(
    bacterias::Bacterias,
    phages::Phages,
    t_max::Int,
    division_rate::Float64,
    death_rate::Float64,
    spacer_loss_rate::Float64,
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
    iter = 0

    snapshots = Vector{Tuple{Float64, Int, Int, Int, Int}}()

    cache = init_rates(bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate, spacer_loss_rate)

    while t < t_max

        if cache.lambda == 0
            break
        end

        num_bacteria = sum(values(bacterias))
        num_phage = sum(values(phages))

        if min(num_bacteria, num_phage) < θ

            tau = rand(Exponential(1.0 / cache.lambda))
            (event_type, event_data) = sample_event(cache)

            if event_type == :division
                apply_division!(bacterias, event_data)
                cache.N_B += 1
                update_rates!(cache, bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate, spacer_loss_rate, :division, event_data)

            elseif event_type == :death
                apply_death!(bacterias, event_data)
                cache.N_B -= 1
                update_rates!(cache, bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate, spacer_loss_rate, :death, event_data)

            elseif event_type == :spacer_loss
                spacers = event_data
                old_bac_keys = Set(keys(bacterias))
                apply_spacer_loss!(bacterias, spacers)
                update_rates!(cache, bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate, spacer_loss_rate, :spacer_loss, spacers)
                for new_spacers in setdiff(keys(bacterias), old_bac_keys)
                    update_rates!(cache, bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate, spacer_loss_rate, :new_clone, new_spacers)
                end

            elseif event_type == :phage_decay
                apply_phage_decay!(phages, event_data)
                update_rates!(cache, bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate, spacer_loss_rate, :phage_decay, event_data)

            elseif event_type == :infection_failed
                spacers, phage_id = event_data
                apply_infection_failed!(bacterias, phages, spacers, phage_id, new_spacer_chance)
                update_rates!(cache, bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate, spacer_loss_rate, :phage_decay, phage_id)

            elseif event_type == :infection_succeeded
                spacers, phage_id = event_data
                old_bac_keys = Set(keys(bacterias))
                old_phage_keys = Set(keys(phages))

                apply_infection_succeeded!(bacterias, phages, spacers, phage_id, burst_size, mutation_chance, new_spacer_chance)

                cache.N_B -= 1

                update_rates!(cache, bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate, spacer_loss_rate, :death, spacers)
                update_rates!(cache, bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate, spacer_loss_rate, :phage_decay, phage_id)

                for new_spacers in setdiff(keys(bacterias), old_bac_keys)
                    cache.N_B += 1
                    update_rates!(cache, bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate, spacer_loss_rate, :new_clone, new_spacers)
                end
                for new_phage_id in setdiff(keys(phages), old_phage_keys)
                    update_rates!(cache, bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate, spacer_loss_rate, :new_phage, new_phage_id)
                end
            end

        else

            tau_calc = (ε_loop * (num_bacteria + num_phage)) / cache.lambda
            tau_max = ε_loop / division_rate
            tau = min(tau_calc, tau_max)

            modified_clones = Set{Set{Int}}()
            modified_phages = Set{Int}()
            new_clones = Set{Set{Int}}()
            new_phages = Set{Int}()

            for (bac_id, rate) in cache.division
                k = rand(Poisson(rate * tau))
                for _ in 1:k
                    if !haskey(bacterias, bac_id) || bacterias[bac_id] <= 0
                        break
                    end
                    apply_division!(bacterias, bac_id)
                    push!(modified_clones, bac_id)
                end
            end

            for (bac_id, rate) in cache.death
                k = rand(Poisson(rate * tau))
                for _ in 1:k
                    if !haskey(bacterias, bac_id) || bacterias[bac_id] <= 0
                        break
                    end
                    apply_death!(bacterias, bac_id)
                    push!(modified_clones, bac_id)
                end
            end

            for (ph_id, rate) in cache.phage_decay
                k = rand(Poisson(rate * tau))
                for _ in 1:k
                    if !haskey(phages, ph_id) || phages[ph_id] <= 0
                        break
                    end
                    apply_phage_decay!(phages, ph_id)
                    push!(modified_phages, ph_id)
                end
            end

            for ((bac_id, phage_id), rate) in cache.infection_failed
                k = rand(Poisson(rate * tau))
                for _ in 1:k
                    if !haskey(bacterias, bac_id) || bacterias[bac_id] <= 0 ||
                       !haskey(phages, phage_id) || phages[phage_id] <= 0
                        break
                    end
                    apply_infection_failed!(bacterias, phages, bac_id, phage_id, new_spacer_chance)
                    push!(modified_phages, phage_id)
                end
            end

            for ((bac_id, phage_id), rate) in cache.infection_succeeded
                k = rand(Poisson(rate * tau))
                old_bac_keys = Set(keys(bacterias))
                old_phage_keys = Set(keys(phages))
                for _ in 1:k
                    if !haskey(bacterias, bac_id) || bacterias[bac_id] <= 0 ||
                       !haskey(phages, phage_id) || phages[phage_id] <= 0
                        break
                    end
                    apply_infection_succeeded!(bacterias, phages, bac_id, phage_id, burst_size, mutation_chance, new_spacer_chance)
                    for new_spacers in setdiff(keys(bacterias), old_bac_keys)
                        push!(new_clones, new_spacers)
                    end
                    for new_phage_id in setdiff(keys(phages), old_phage_keys)
                        push!(new_phages, new_phage_id)
                    end
                    push!(modified_clones, bac_id)
                    push!(modified_phages, phage_id)
                    old_bac_keys = Set(keys(bacterias))
                    old_phage_keys = Set(keys(phages))
                end
            end

            for (spacers, rate) in cache.spacer_loss
                k = rand(Poisson(rate * tau))
                old_bac_keys = Set(keys(bacterias))
                for _ in 1:k
                    if !haskey(bacterias, spacers) || bacterias[spacers] <= 0 || isempty(spacers)
                        break
                    end
                    apply_spacer_loss!(bacterias, spacers)
                    push!(modified_clones, spacers)
                end
                for new_spacers in setdiff(keys(bacterias), old_bac_keys)
                    push!(new_clones, new_spacers)
                end
            end

            cache.N_B = sum(values(bacterias))

            for spacers in new_clones
                update_rates!(cache, bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate, spacer_loss_rate, :new_clone, spacers)
            end
            for phage_id in new_phages
                update_rates!(cache, bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate, spacer_loss_rate, :new_phage, phage_id)
            end
            for spacers in modified_clones
                update_rates!(cache, bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate, spacer_loss_rate, :division, spacers)
            end
            for phage_id in modified_phages
                update_rates!(cache, bacterias, phages, division_rate, death_rate, K, phage_decay, infection_rate, spacer_loss_rate, :phage_decay, phage_id)
            end

        end

        t += tau

        iter += 1
        if iter % 1000 == 0
            cache.lambda = sum(values(cache.division)) +
                           sum(values(cache.death)) +
                           sum(values(cache.spacer_loss)) +
                           sum(values(cache.phage_decay)) +
                           sum(values(cache.infection_failed)) +
                           sum(values(cache.infection_succeeded))
            cache.N_B = sum(values(bacterias))
        end

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