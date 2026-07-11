using Test
include("../src/CHASE.jl")
using .CHASE

@testset "recognize" begin
    @test recognize(Set([1, 2, 3]), 1) == true
    @test recognize(Set([1, 2, 3]), 5) == false
    @test recognize(Set{Int}(), 1)     == false
end

@testset "apply_division!" begin
    bacterias = Bacterias(Set([1, 2]) => 10)
    apply_division!(bacterias, Set([1, 2]))
    @test bacterias[Set([1, 2])] == 11
end

@testset "apply_death!" begin
    bacterias = Bacterias(Set([1, 2]) => 3)
    apply_death!(bacterias, Set([1, 2]))
    @test bacterias[Set([1, 2])] == 2

    bacterias2 = Bacterias(Set([1]) => 1)
    apply_death!(bacterias2, Set([1]))
    @test !haskey(bacterias2, Set([1]))
end

@testset "apply_phage_decay!" begin
    phages = Phages(42 => 3)
    apply_phage_decay!(phages, 42)
    @test phages[42] == 2

    phages2 = Phages(42 => 1)
    apply_phage_decay!(phages2, 42)
    @test !haskey(phages2, 42)
end

@testset "apply_infection_succeeded!" begin
    bacterias = Bacterias(Set([1]) => 5)
    phages = Phages(42 => 10)
    apply_infection_succeeded!(bacterias, phages, Set([1]), 42, 100, 0.0, 0.0)
    @test bacterias[Set([1])] == 4
    @test phages[42] == 110
    @test length(phages) == 1

    bacterias2 = Bacterias(Set([1]) => 5)
    phages2 = Phages(42 => 10)
    apply_infection_succeeded!(bacterias2, phages2, Set([1]), 42, 100, 0.0, 1.0)
    @test haskey(bacterias2, Set([1, 42]))
    @test bacterias2[Set([1, 42])] == 1
    @test bacterias2[Set([1])] == 5
    @test phages2[42] == 9
end

@testset "apply_infection_failed!" begin
    bacterias = Bacterias(Set([1]) => 5)
    phages = Phages(42 => 1)
    apply_infection_failed!(bacterias, phages, Set([1]), 42, 0.0)
    @test !haskey(phages, 42)
    @test bacterias[Set([1])] == 5

    bacterias2 = Bacterias(Set([1]) => 5)
    phages2 = Phages(42 => 3)
    apply_infection_failed!(bacterias2, phages2, Set([1]), 42, 1.0)
    @test phages2[42] == 2
end

@testset "init_rates" begin
    bacterias = Bacterias(Set([1]) => 100)
    phages = Phages(42 => 20)
    cache = init_rates(bacterias, phages, 0.1, 0.05, 10000, 0.001, 0.001)
    @test cache.lambda > 0.0
    @test haskey(cache.division, Set([1]))
    @test haskey(cache.death, Set([1]))
    @test haskey(cache.phage_decay, 42)
    @test haskey(cache.infection_succeeded, (Set([1]), 42))

    bacterias2 = Bacterias()
    phages2 = Phages()
    cache2 = init_rates(bacterias2, phages2, 0.1, 0.05, 10000, 0.001, 0.001)
    @test cache2.lambda == 0.0
end

@testset "gillespie exact" begin
    bacterias = Bacterias(Set{Int}() => 10)
    phages = Phages(1 => 5)
    snapshots = gillespie(bacterias, phages, 10, 0.3, 0.1, 0.001, 10000, 50, 0.0001, 0.00001, 0.01, 0, 0.03)
    @test length(snapshots) > 0
    @test length(snapshots[1]) == 5
    @test snapshots[1][1] isa Float64
    @test snapshots[1][2] isa Int
end

@testset "gillespie tau-leaping" begin
    bacterias = Bacterias(Set{Int}() => 10)
    phages = Phages(1 => 5)
    snapshots = gillespie(bacterias, phages, 10, 0.3, 0.1, 0.001, 10000, 50, 0.0001, 0.00001, 0.01, 1_000_000, 0.03)
    @test length(snapshots) > 0
    @test length(snapshots[1]) == 5
    @test snapshots[1][1] isa Float64
    @test snapshots[1][2] isa Int
end