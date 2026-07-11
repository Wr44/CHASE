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

@testset "apply_infection_succeeded!" begin
    bacterias = Bacterias(Set([1]) => 5)
    phages = Phages(42 => 10)
    apply_infection_succeeded!(bacterias, phages, Set([1]), 42, 100, 0.0)
    @test bacterias[Set([1])] == 4
    @test phages[42] == 110
    @test length(phages) == 1
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
    @test haskey(bacterias2, Set([1, 42]))
    @test bacterias2[Set([1, 42])] == 1
    @test bacterias2[Set([1])] == 4
end

@testset "calculate_rates" begin
    bacterias = Bacterias(Set([1]) => 100)
    phages = Phages(42 => 20)
    (rates, lambda) = calculate_rates(bacterias, phages, 0.1, 0.05, 0.001)
    @test lambda > 0.0
    @test length(rates) == 3

    bacterias2 = Bacterias()
    phages2 = Phages()
    (rates2, lambda2) = calculate_rates(bacterias2, phages2, 0.1, 0.05, 0.001)
    @test lambda2 == 0.0
end