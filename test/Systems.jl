@testset "Equivalenz TwoLevelChain - DissipativeLambaChain" begin
	# ----------------------------------------
	# At zero γd both models must coincide
	Γ = rand()

  problemA = WaveguideProblem(TwoLevelChain(1, Γ = Γ),
		WavePacket(SoftBoxMode(τ = 1.0), Coherent(1.0)), 2.
	)
  solA     = solve(problemA, reltol = 1e-10)[2]

  problemB = WaveguideProblem(DissipativeLambdaChain(1, Γ = Γ, γd = 0.0),
		WavePacket(SoftBoxMode(τ = 1.0), Coherent(1.0)), 2.
	)
  solB     = solve(problemB, reltol = 1e-10)[2]

  # Project on the groundstate
  σ_GG_A = transition(NLevelBasis(2), 1, 1)
  σ_GG_B = transition(NLevelBasis(3), 1, 1)
  @test sum(abs, expect(σ_GG_A, solA) .- expect(σ_GG_B, solB) ) < 1e-5

  # Project on the excited state
  σ_WW_A = transition(NLevelBasis(2), 2, 2)
  σ_WW_B = transition(NLevelBasis(3), 2, 2)
  @test sum(abs, expect(σ_WW_A, solA) .- expect(σ_WW_B, solB) ) < 1e-5
end

@testset "TwoLevelAtoms - Chirality" begin
	# ----------------------------------------
	# The dynamics of an atom is not impacted by latter atoms
  Γ = rand()
  σ_WW_i = transition(NLevelBasis(2), 2, 2)

  problem1 = WaveguideProblem(
    TwoLevelChain(1, Γ = Γ), ContinuousWave(1.0), 3.0
  )
  sol1 = solve(problem1, reltol = 1e-10)[2]

  problem2 = WaveguideProblem(
    TwoLevelChain(2, Γ = Γ), ContinuousWave(1.0), 3.0
  )
  sol2 = solve(problem2, reltol = 1e-10)[2]
	sol21, sol22 = ptrace.(sol2, 2), ptrace.(sol2, 1)

  problem3 = WaveguideProblem(
    TwoLevelChain(3, Γ = Γ), ContinuousWave(1.0), 3.0
  )
  sol3 = solve(problem3, reltol = 1e-10)[2]
	sol31, sol32 = [ptrace(sol, [2, 3]) for sol ∈ sol3], [ptrace(sol, [1, 3]) for sol ∈ sol3]

	@test sum(abs, tracedistance.( sol1,  sol21 )) < 1e-5
	@test sum(abs, tracedistance.( sol1,  sol31 )) < 1e-5
	@test sum(abs, tracedistance.( sol22, sol32 )) < 1e-5
end

@testset "DissipativeLambdaChain - Chirality" begin
	# ----------------------------------------
	# The dynamics of an atom is not impacted by latter atoms
  Γ, γd = rand(), rand()
  σ_WW_i = transition(NLevelBasis(3), 2, 2)

  problem1 = WaveguideProblem(
    DissipativeLambdaChain(1, γd = γd, Γ = Γ), ContinuousWave(1.0), 3.0
  )
  sol1 = solve(problem1, reltol = 1e-10)[2]

  problem2 = WaveguideProblem(
    DissipativeLambdaChain(2, γd = γd, Γ = Γ), ContinuousWave(1.0), 3.0
  )
  sol2 = solve(problem2, reltol = 1e-10)[2]
	sol21, sol22 = ptrace.(sol2, 2), ptrace.(sol2, 1)

  problem3 = WaveguideProblem(
    DissipativeLambdaChain(3, γd = γd, Γ = Γ), ContinuousWave(1.0), 3.0
  )
  sol3 = solve(problem3, reltol = 1e-10)[2]
	sol31, sol32 = [ptrace(sol, [2, 3]) for sol ∈ sol3], [ptrace(sol, [1, 3]) for sol ∈ sol3]

  @test sum(abs, tracedistance.( sol1,  sol21 )) < 1e-5
	@test sum(abs, tracedistance.( sol1,  sol31 )) < 1e-5
	@test sum(abs, tracedistance.( sol22, sol32 )) < 1e-5
end
