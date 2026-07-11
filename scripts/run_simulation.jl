include("../src/CHASE.jl")
using .CHASE

β = 0.1
δ = 0.05
φ = 0.001
B = 100
μ = 0.01
ε = 0.1

tmax = 100

bacterias = Bacterias(Set{Int}() => 100)
phages = Phages(1 => 20)

snapshots = gillespie(
    bacterias,
    phages,
    tmax,
    β, δ, B, φ, μ, ε
)

for (t, bact, p) in snapshots 
    println("t = $t | bacterias = $bact | phages = $p")
end 