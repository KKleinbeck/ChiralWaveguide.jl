using SpecialFunctions
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
    Mode(x -> twoAtomOut(-x); compression = Exponential(σ = σ), reltol = 1e-10),
    (-5σ, 20σ)
  )
  sol = ptrace(solve(problem, reltol = 1e-10)[2][end], [1, 2, 3])

  # Transfer chance
  n̂ = number(basis(sol))
  @test real( 1 - expect(n̂, sol) ) < 1e-5
end

@testset "Fock State scattering - Kiilerich" begin
	# ----------------------------------------
	# Recalculate Figure 2 of PRL 123, 123604 (2019)
	function outputModeFunctionFig2(t, γ = 1.0, τ = 4.0)
		abs(t-τ) > 10 && return 0.0
	  return exp(-(t - τ)^2/2) -
	    γ * √(π/2) * exp(-γ * (t - τ) / 2 + γ^2 / 8) * erfc((-2(t - τ) + γ) / (2 * √(2)))
	end

	WP = WavePacket(GaussMode(τ = 4.0), Fock(1))
	outputMode = Mode(t -> outputModeFunctionFig2(t, 1.0, 4.0), compression = :algebraic)

	problemA = WaveguideProblem(TwoLevelChain(1), WP, outputMode, 13.0)
	solA     = solve(problemA, reltol = 1e-10)[2]
	b        = basis(solA[end])

	n̂ᵢ, n̂ₒ = embed(b, 1, number(b.bases[1])), embed(b, 3, number(b.bases[3]))
	σ⁺σ⁻   = embed(b, 2, transition(NLevelBasis(2), 2, 2))
	excitations = n̂ᵢ + n̂ₒ + σ⁺σ⁻

	# particle Number conservation and correct final state
	particleNumbers = expect(excitations, solA) .|> real
	@test isapprox.(particleNumbers, 1.0, atol = 1e-5) |> all
	@test isapprox(ptrace(solA[end], [1, 2]).data[end, end], 1.0, atol = 1e-3)

	# ----------------------------------------
	# Two photon scattering
	WP2 = WavePacket(GaussMode(τ = 4.0), Fock(2))

	problemB = WaveguideProblem(TwoLevelChain(1), WP2, WP2.mode, 13.0)
	solB     = solve(problemB, reltol = 1e-10)[2]

	@test isapprox(
		ptrace(solB[end], [1, 2]).data[end, end] |> real,
		0.639387^2, atol = 1e-3
	) # See `Two Photon Overlaps.nb` for magic numbers

	problemB = WaveguideProblem(TwoLevelChain(1), WP2, outputMode, 13.0)
	solB     = solve(problemB, reltol = 1e-10)[2]

	@test isapprox(
		ptrace(solB[end], [1, 2]).data[end, end] |> real,
		0.17789^2, atol = 1e-3
	) # See `Two Photon Overlaps.nb` for magic numbers
end

@testset "Perfect Absorption" begin
	# ----------------------------------------
	# With a suitable state perfect absorption is almost possible
  γd =  0.5
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

@testset "Exact Coherent Evolution" begin
	# ----------------------------------------
	# Compare against the exact time-evolution of coherently driven two level atom

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
