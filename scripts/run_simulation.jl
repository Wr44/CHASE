include("../src/CHASE.jl")

using Plots
using .CHASE

β = 1.0
δ = 0.1
φ = 1e-7
B = 100
μ = 1e-5
ε = 1e-4
θ = 15
ε_loop = 0.03 

tmax = 500

bacterias = Bacterias(Set{Int}() => 1000)
phages = Phages(1 => 10000)

snapshots = gillespie(
    bacterias,
    phages,
    tmax,
    β, δ, B, φ, μ, ε, θ, ε_loop
)

t_vals = [s[1] for s in snapshots]
bact_vals = [s[2] for s in snapshots]
phage_vals = [s[3] for s in snapshots]
clones = [s[4] for s in snapshots]
strains = [s[5] for s in snapshots]

p1 = plot(t_vals, bact_vals, label="Bacteria", ylabel="Count")
plot!(twinx(), t_vals, phage_vals, label="Phages", color=:red, ylabel="Phages")

p2 = plot(t_vals, clones, label="Clones bactériens", ylabel="Diversité", color=:blue)
plot!(t_vals, strains, label="Souches phagiques", color=:red)

plot(p1, p2, layout=(2,1), xlabel="Time", title="CHASE - Coevolution Dynamics")
savefig("recap.png")