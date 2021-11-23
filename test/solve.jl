@testset "Specifying Initial state - Continuous Wave" begin
	# ----------------------------------------
	# Manually specifying the initial state for ContinuousWave
  problem = WaveguideProblem(TwoLevelChain(1), ContinuousWave(1.0), 10.0)
  solA = solve(problem, reltol = 1e-10)[2]

  # Project on the Rydberg manifold
  σ_Ryd = transition(NLevelBasis(2), 2, 2)

	ρ₀ = nlevelstate(NLevelBasis(2), 1)
  solB = solve(problem, ρ₀ = ρ₀, reltol = 1e-10)[2]

  @test sum( tracedistance.(solA, solB) ) < 1e-4
end

@testset "Specifying Initial state - Fock" begin
	# ----------------------------------------
	# Variant of "Systems/Perfect Absorption" test
  problem = WaveguideProblem(TwoLevelChain(1), WavePacket(ExpMode(), Fock(1)), (-10.0, 0.0))
  solA = solve(problem, reltol = 1e-10)[2]

  # Project on the Rydberg manifold
  σ_Ryd = transition(NLevelBasis(2), 2, 2)

	ρ₀ = fockstate(FockBasis(1), 1) ⊗ nlevelstate(NLevelBasis(2), 1)
  solB = solve(problem, ρ₀ = ρ₀, reltol = 1e-10)[2]

  @test sum( tracedistance.(solA, solB) ) < 1e-4
end

@testset "Overwriting Nouts" begin
	problem = WaveguideProblem(TwoLevelChain(1), WavePacket(ExpMode(), Fock(1)), ExpMode(), 10.0)
	basisF, basisA = FockBasis(1), NLevelBasis(2)
	ρ₀ = fockstate(basisF, 1) ⊗ nlevelstate(basisA, 1) ⊗ fockstate(basisF, 0)
	@test_logs (:warn, "`Nouts` and `ρ₀` simultaneously set. Ignoring `Nouts`") begin
		solve(problem, ρ₀ = ρ₀, Nouts = [3])
		true
	end
end
