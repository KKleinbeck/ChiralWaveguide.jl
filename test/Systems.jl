using SpecialFunctions

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

@testset "Perfect Absorption" begin
	# ----------------------------------------
	# With a suitable state perfect absorption is almost possible
  γd =  0.65*rand() + 0.3 # smaller γd need more careful numerical treatment
	γd =  0.3
  t₀ = -10.0/(1 - γd)

  problemA = WaveguideProblem(
		DissipativeLambdaChain(1, γd = γd),
    WavePacket(ExpMode(γ = 1 - γd), Fock(1)),
    (t₀, 0.0)
  )
  solA = ptrace.(solve(problemA, reltol = 1e-10)[2], 1)[end]

  # Project on the Rydberg manifold
  σ_Ryd = transition(NLevelBasis(3), 2, 2) + transition(NLevelBasis(3), 3, 3)
  @test real( 1 - expect(σ_Ryd, solA) ) < 1e-4


	problemB = WaveguideProblem(
		DissipativeLambdaChain(1, γd = γd),
		WavePacket(ExpMode(γ = 1 - γd), Fock(1)),
		(t₀, 10.0)
	)
  solB = ptrace.(solve(problemB, reltol = 1e-10)[2], 1)[end]

  # Project on the groundstate
  σ_GG = transition(NLevelBasis(3), 1, 1)
  @test abs( (1 - γd) / (1 + γd) - expect(σ_GG, solB) ) < 1e-4

  # Project on the dark state
  σ_DD = transition(NLevelBasis(3), 3, 3)
  @test abs( 2γd / (1 + γd) - expect(σ_DD, solB) ) < 1e-4
end

@testset "Two Atom scattering" begin
	# ----------------------------------------
	# Exact result for scattering a Gaussian mode
  σ = 1 + rand() # smaller σ need more careful numerical treatment
  function twoAtomOut(x)
    norm = √(√(π)*σ)

    return (1 + σ^2) * exp(-x^2/(2σ^2)) / norm -
      (norm / √(2)) * (x + 2 + σ^2 / 2) * exp((4x + σ^2)/8) * erfc((2x + σ^2) / (2 * √(2) * σ))
  end

  problem = WaveguideProblem(
    TwoLevelChain(2),
    WavePacket(GaussMode(σ = σ), Fock(1)),
    Mode(x -> twoAtomOut(-x); compression = Exponential(σ = σ)),
    (-5σ, 20σ)
  )
  sol = ptrace(solve(problem, reltol = 1e-10)[2][end], [1, 2, 3])

  # Transfer chance
  n̂ = number(sol.basis_l)
  @test real( 1 - expect(n̂, sol) ) < 1e-5
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

@testset "Exact Coherent Evolution" begin
	# ----------------------------------------
	# Compare against the exact time-evolution

	function blochVectorToDensity(vec)
	  ρ = nlevelstate(NLevelBasis(2), 1) |> dm

	  ρ.data[1, 1] = vec[1]
	  ρ.data[1, 2] = vec[2] / 2
	  ρ.data[2, 1] = conj(vec[2]) / 2
	  ρ.data[2, 2] = vec[3]

	  return ρ
	end

	function exactDynamics(t, κ, α; transform = false)
	  Ω = 4 * √(κ) * α
	  ω = √(4*Ω^2 - κ^2)

	  ρSS = [Ω^2 + 4κ^2, 4im * κ * Ω, Ω^2] ./ (2Ω^2 + 4κ^2)
	  ρ1  = [-1, (1.0im * κ + ω) / Ω, 1]
	  ρ2  = [-1, (1.0im * κ - ω) / Ω, 1]

	  c1 = Ω^2 / (-ω^2 + 3im * ω * κ) * exp((-3κ - 1.0im * ω) * t / 4)
	  c2 = Ω^2 / (-ω^2 - 3im * ω * κ) * exp((-3κ + 1.0im * ω) * t / 4)

	  ρ = ρSS .+ c1 .* ρ1 .+ c2 .* ρ2 |> blochVectorToDensity

	  if transform
	    ρ.data[1,2] *=  1im
	    ρ.data[2,1] *= -1im
	  end
	  return ρ
	end

	# ----------------------------------------
	# Test time evolution and steady state of superatom
  id, σWG, σGW = one(NLevelBasis(2)), transition(NLevelBasis(2), 2, 1),
		transition(NLevelBasis(2), 1, 2)
	αs = range(0.2, 2.0, length = 19)
	for αi ∈ αs
		problem = WaveguideProblem(TwoLevelChain(1), ContinuousWave(αi), 20.0)
		ts, sol = solve(problem, reltol = 1e-10)

		ρAtomsExact = exactDynamics.(ts, 1.0, αi, transform = true)
		@test sum(abs, tracedistance.(sol, ρAtomsExact)) / length(ts) < 1e-7

		@test abs(expect((αi*id + σWG) * (αi*id + σGW), sol[end]) - αi^2) < 1e-6
	end
end
