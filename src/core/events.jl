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

function apply_infection_succeeded!(
    bacterias::Bacterias, 
    phages::Phages, 
    spacers::Set{Int}, 
    phage_id::Int, 
    burst_size::Int, 
    mutation_chance::Float64
    )
    apply_death!(bacterias, spacers)

    d = Binomial(burst_size, mutation_chance) 
    n = rand(d)

    phages[phage_id] += burst_size - n 

    for _ in 1:n 

        new_id = rand(Int) 

        while haskey(phages, new_id)
            new_id = rand(Int)
        end

        phages[new_id] = 1 
    end 
end

function apply_infection_failed!(
    bacterias::Bacterias, 
    phages::Phages, 
    spacers::Set{Int}, 
    phage_id::Int, 
    new_spacer_chance::Float64
    )
    
    if (phages[phage_id] == 1) 
        delete!(phages, phage_id)
    else 
        phages[phage_id] -= 1 
    end 

    if (rand() < new_spacer_chance)
        new_spacers = union(spacers, Set([phage_id]))
        apply_death!(bacterias, spacers)

        if haskey(bacterias, new_spacers)
             bacterias[new_spacers] += 1
        else
            bacterias[new_spacers] = 1 
        end
    end

end
