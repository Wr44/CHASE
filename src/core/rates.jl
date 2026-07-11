mutable struct RatesCache
    division :: Dict{Set{Int}, Float64}
    death :: Dict{Set{Int}, Float64}
    phage_decay :: Dict{Int, Float64}
    infection_failed :: Dict{Tuple{Set{Int}, Int}, Float64}
    infection_succeeded :: Dict{Tuple{Set{Int}, Int}, Float64}
    lambda :: Float64
end 

function init_rates(
    bacterias::Bacterias,
    phages::Phages,
    division_rate::Float64,
    death_rate::Float64,
    K::Int,
    phage_decay::Float64,
    infection_rate::Float64,
)::RatesCache

    rates = RatesCache(
        Dict(),Dict(),Dict(),Dict(),Dict(),0.0
    )

    N_B = sum(values(bacterias)) 
    
    for (spacers, bacteria_count) in bacterias

        rates.division[spacers] = division_rate * bacteria_count
        rates.death[spacers] = (death_rate + N_B/K)* bacteria_count

        for (phage_id, phage_count) in phages
            if recognize(spacers, phage_id)
                rates.infection_failed[(spacers, phage_id)] = infection_rate*bacteria_count*phage_count
            else
                rates.infection_succeeded[(spacers, phage_id)] = infection_rate*bacteria_count*phage_count
            end    
        end
    end

    for (phage_id, phage_count) in phages
        rates.phage_decay[phage_id] = phage_decay * phage_count
    end

    rates.lambda = sum(values(rates.division)) +
               sum(values(rates.death)) +
               sum(values(rates.phage_decay)) +
               sum(values(rates.infection_failed)) +
               sum(values(rates.infection_succeeded))
    return rates
end


function update_rates!(
    cache::RatesCache,
    bacterias::Bacterias,
    phages::Phages,
    division_rate::Float64,
    death_rate::Float64,
    K::Int,
    phage_decay_rate::Float64,
    infection_rate::Float64,
    event_type::Symbol,
    event_data::Any
)
    N_B = sum(values(bacterias))

    if event_type == :division || event_type == :death
        spacers = event_data
        n = get(bacterias, spacers, 0)

        if n == 0
            old_div = get(cache.division, spacers, 0.0)
            old_dth = get(cache.death, spacers, 0.0)
            cache.lambda -= old_div + old_dth
            delete!(cache.division, spacers)
            delete!(cache.death, spacers)
            for (phage_id, _) in phages
                key = (spacers, phage_id)
                cache.lambda -= get(cache.infection_failed, key, 0.0) +
                                get(cache.infection_succeeded, key, 0.0)
                delete!(cache.infection_failed, key)
                delete!(cache.infection_succeeded, key)
            end
        else
            old_div = get(cache.division, spacers, 0.0)
            old_dth = get(cache.death, spacers, 0.0)
            new_div = division_rate * n
            new_dth = (death_rate + N_B / K) * n
            cache.division[spacers] = new_div
            cache.death[spacers] = new_dth
            cache.lambda += (new_div - old_div) + (new_dth - old_dth)

            for (phage_id, phage_count) in phages
                key = (spacers, phage_id)
                old_inf = get(cache.infection_failed, key, 0.0) +
                          get(cache.infection_succeeded, key, 0.0)
                new_inf = infection_rate * n * phage_count
                if recognize(spacers, phage_id)
                    cache.infection_failed[key] = new_inf
                    delete!(cache.infection_succeeded, key)
                else
                    cache.infection_succeeded[key] = new_inf
                    delete!(cache.infection_failed, key)
                end
                cache.lambda += new_inf - old_inf
            end
        end

    elseif event_type == :phage_decay
        phage_id = event_data
        m = get(phages, phage_id, 0)

        if m == 0
            old_decay = get(cache.phage_decay, phage_id, 0.0)
            cache.lambda -= old_decay
            delete!(cache.phage_decay, phage_id)
            for (spacers, _) in bacterias
                key = (spacers, phage_id)
                cache.lambda -= get(cache.infection_failed, key, 0.0) +
                                get(cache.infection_succeeded, key, 0.0)
                delete!(cache.infection_failed, key)
                delete!(cache.infection_succeeded, key)
            end
        else
            old_decay = get(cache.phage_decay, phage_id, 0.0)
            new_decay = phage_decay_rate * m
            cache.phage_decay[phage_id] = new_decay
            cache.lambda += new_decay - old_decay

            for (spacers, bacteria_count) in bacterias
                key = (spacers, phage_id)
                old_inf = get(cache.infection_failed, key, 0.0) +
                          get(cache.infection_succeeded, key, 0.0)
                new_inf = infection_rate * bacteria_count * m
                if recognize(spacers, phage_id)
                    cache.infection_failed[key] = new_inf
                    delete!(cache.infection_succeeded, key)
                else
                    cache.infection_succeeded[key] = new_inf
                    delete!(cache.infection_failed, key)
                end
                cache.lambda += new_inf - old_inf
            end
        end

    elseif event_type == :infection_succeeded || event_type == :infection_failed
        spacers, phage_id = event_data
        update_rates!(cache, bacterias, phages, division_rate, death_rate, K,
                      phage_decay_rate, infection_rate, :division, spacers)
        update_rates!(cache, bacterias, phages, division_rate, death_rate, K,
                      phage_decay_rate, infection_rate, :phage_decay, phage_id)

    elseif event_type == :new_clone
        spacers = event_data
        n = get(bacterias, spacers, 0)
        new_div = division_rate * n
        new_dth = (death_rate + N_B / K) * n
        cache.division[spacers] = new_div
        cache.death[spacers] = new_dth
        cache.lambda += new_div + new_dth
        for (phage_id, phage_count) in phages
            key = (spacers, phage_id)
            new_inf = infection_rate * n * phage_count
            if recognize(spacers, phage_id)
                cache.infection_failed[key] = new_inf
            else
                cache.infection_succeeded[key] = new_inf
            end
            cache.lambda += new_inf
        end

    elseif event_type == :new_phage
        phage_id = event_data
        m = get(phages, phage_id, 0)
        new_decay = phage_decay_rate * m
        cache.phage_decay[phage_id] = new_decay
        cache.lambda += new_decay
        for (spacers, bacteria_count) in bacterias
            key = (spacers, phage_id)
            new_inf = infection_rate * bacteria_count * m
            if recognize(spacers, phage_id)
                cache.infection_failed[key] = new_inf
            else
                cache.infection_succeeded[key] = new_inf
            end
            cache.lambda += new_inf
        end
    end
end