using Distributions

function apply_division!(
    bacterias::Bacterias, 
    spacers::Set{Int}
    )
    bacterias[spacers] += 1 
end 


function apply_death!(
    bacterias::Bacterias, 
    spacers::Set{Int}
    )
    if (bacterias[spacers] == 1) 
        delete!(bacterias, spacers)
    else 
        bacterias[spacers] -= 1 
    end 
end 

function apply_phage_decay!(
    phages::Phages,
    phage_id::Int
)
    if (phages[phage_id] == 1) 
        delete!(phages, phage_id)
    else 
        phages[phage_id] -= 1 
    end 
end

function apply_infection_succeeded!(
    bacterias::Bacterias,
    phages::Phages,
    spacers::Set{Int},
    phage_id::Int,
    burst_size::Int,
    mutation_chance::Float64,
    new_spacer_chance::Float64
)
    if rand() < new_spacer_chance
        new_spacers = union(spacers, Set([phage_id]))
        if haskey(bacterias, new_spacers)
            bacterias[new_spacers] += 1
        else
            bacterias[new_spacers] = 1
        end
        apply_phage_decay!(phages, phage_id)
    else
        apply_death!(bacterias, spacers)
        d = Binomial(burst_size, mutation_chance)
        n = rand(d)
        phages[phage_id] = get(phages, phage_id, 0) + burst_size - n
        for _ in 1:n
            new_id = rand(Int)
            while haskey(phages, new_id)
                new_id = rand(Int)
            end
            phages[new_id] = 1
        end
    end
end

function apply_infection_failed!(
    phages::Phages,
    phage_id::Int,
)
    apply_phage_decay!(phages, phage_id)
end


function apply_spacer_loss!(
    bacterias::Bacterias,
    spacers::Set{Int}
) 
    new_spacers = copy(spacers)
    removed = rand(spacers)

    delete!(new_spacers, removed)

    if haskey(bacterias, new_spacers) 
        bacterias[new_spacers] += 1
    else
        bacterias[new_spacers] = 1
    end

    if (bacterias[spacers] == 1) 
        delete!(bacterias, spacers)
    else 
        bacterias[spacers] -= 1 
    end 

end
