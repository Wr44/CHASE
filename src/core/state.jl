const Bacterias = Dict{Set{Int}, Int}
const Phages = Dict{Int, Int}

function recognize(
    spacers::Set{Int},
    phage_id::Int
):: Bool
    return phage_id in spacers
end

