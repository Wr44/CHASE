include("../src/CHASE.jl")

using Plots
using .CHASE

β = 1.0
δ = 0.1
K = 50_000
φ = 1e-7
B = 100
d = 0.02
μ = 1e-6
ε = 1e-4
κ = 1e-5
θ = 15
ε_loop = 0.03

tmax = 1000

bacterias = Bacterias(Set{Int}() => 1000)
phages = Phages(1 => 1000)


snapshots = gillespie(
    bacterias,
    phages,
    tmax,
    β, δ, κ, d, K, B, φ, μ, ε, θ, ε_loop
)

step_size = max(1, length(snapshots) ÷ 7000) 
snapshots_light = snapshots[1:step_size:end]

t_vals = [s[1] for s in snapshots_light]
bact_vals = [s[2] for s in snapshots_light]
phage_vals = [s[3] for s in snapshots_light]
clones = [s[4] for s in snapshots_light]
strains = [s[5] for s in snapshots_light]

p1 = plot(t_vals, bact_vals, label="Bacteria", ylabel="Count")
plot!(twinx(), t_vals, phage_vals, label="Phages", color=:red, ylabel="Phages")

p2 = plot(t_vals, clones, label="Spacers", ylabel="Spacers", color=:blue)
plot!(twinx(), t_vals, strains, label="Phages types", color=:red, ylabel="Phages types")

plot(p1, p2, layout=(2,1), xlabel="Time", title="CHASE - Coevolution Dynamics")
savefig("recap.png")