@testset "DrivenProblem - Absence of Dissipator" begin
	# ----------------------------------------
	# DrivenProblem - empty Dissipator
  testBasis = SpinBasis(1//2)
  σm, σp = sigmam(testBasis), sigmap(testBasis)

	WP, mode = WavePacket(SoftBoxMode(), Coherent(1.0)), SoftBoxMode()

  problemA = WaveguideProblem((σm+σp,     σm, spindown(testBasis)), WP, 1.)
  solA     = solve(problemA)[2]

  problemB = WaveguideProblem((σm+σp, [], σm, spindown(testBasis)), WP, 1.)
  solB     = solve(problemB)[2]

  @test sum(tracedistance.(solA, solB)) ≈ 0
end

@testset "ScatterProblem - Absence of Dissipator" begin
	# ----------------------------------------
	# ScatterProblem - empty Dissipator
  testBasis = SpinBasis(1//2)
  σm, σp = sigmam(testBasis), sigmap(testBasis)

	WP, mode = WavePacket(SoftBoxMode(), Coherent(1.0)), SoftBoxMode()

  problemA = WaveguideProblem((σm+σp,     σm, spindown(testBasis)), WP, mode, 1.)
  solA     = solve(problemA)[2]

  problemB = WaveguideProblem((σm+σp, [], σm, spindown(testBasis)), WP, mode, 1.)
  solB     = solve(problemB)[2]

  @test sum(tracedistance.(solA, solB)) ≈ 0
end

@testset "DrivenProblem - ContinuousWave" begin
	# ----------------------------------------
	# ScatterProblem - empty Dissipator
  testBasis = SpinBasis(1//2)
  σm, σp = sigmam(testBasis), sigmap(testBasis)

  problemA = WaveguideProblem((σm+σp, σm, spindown(testBasis)), ContinuousWave(     1.0), 1.)
  solA     = solve(problemA)[2]

  problemB = WaveguideProblem((σm+σp, σm, spindown(testBasis)), ContinuousWave(t -> 1.0), 1.)
  solB     = solve(problemB)[2]

  @test sum(tracedistance.(solA, solB)) ≈ 0
end

@testset "ScatterProblem - ContinuousWave" begin
	# ----------------------------------------
	# ScatterProblem - empty Dissipator
  testBasis = SpinBasis(1//2)
  σm, σp = sigmam(testBasis), sigmap(testBasis)

	mode = SoftBoxMode()

  problemA = WaveguideProblem((σm+σp, σm, spindown(testBasis)), ContinuousWave(     1.0), mode, 1.)
  solA     = solve(problemA)[2]

  problemB = WaveguideProblem((σm+σp, σm, spindown(testBasis)), ContinuousWave(t -> 1.0), mode, 1.)
  solB     = solve(problemB)[2]

  @test sum(tracedistance.(solA, solB)) ≈ 0
end

@testset "Chirality of ScatterProblem" begin
	# ----------------------------------------
	# DrivenProblem - empty Dissipator
  testBasis = SpinBasis(1//2)
  σm, σp = sigmam(testBasis), sigmap(testBasis)

	WP, mode = WavePacket(SoftBoxMode(), Coherent(1.0)), SoftBoxMode()

	problemA = WaveguideProblem((σm+σp, [], σm, spindown(testBasis)), WP, 1.)
  solA     = solve(problemA, reltol = 1e-10)[2]

  problemB = WaveguideProblem((σm+σp, [], σm, spindown(testBasis)), WP, mode, 1.)
  solB     = ptrace.(solve(problemB, reltol = 1e-10)[2], 2)

  @test sum(tracedistance.(solA, solB)) < 1e-5
end

@testset "DrivenProblem - Mollow transform" begin
	# ----------------------------------------
	# Compare the Mollow transform to a manual coherent state
  testBasis = SpinBasis(1//2)
  σm, σp = sigmam(testBasis), sigmap(testBasis)

  problemA = WaveguideProblem((σm+σp, [], σm, spindown(testBasis)),
    WavePacket(SoftBoxMode(τ = 1.0), ArbitraryState(coherentstate(FockBasis(12), 1.0).data)), 3.
  )
  solA     = ptrace.(solve(problemA, reltol = 1e-10)[2], 1)

  problemB = WaveguideProblem((σm+σp, [], σm, spindown(testBasis)),
		WavePacket(SoftBoxMode(τ = 1.0), Coherent(1.0)), 3.
	)
  solB     = solve(problemB, reltol = 1e-10)[2]

  @test sum(tracedistance.(solA, solB)) < 1e-5
end

@testset "ScatterProblem - Mollow transform" begin
	# ----------------------------------------
	# Compare the Mollow transform to a manual coherent state
  testBasis = SpinBasis(1//2)
  σm, σp = sigmam(testBasis), sigmap(testBasis)

  problemA = WaveguideProblem((σm+σp, [], σm, spindown(testBasis)),
    WavePacket(SoftBoxMode(τ = 1.0), ArbitraryState(coherentstate(FockBasis(12), 1.0).data)),
		SoftBoxMode(τ = 2.0), 3.
  )
  solA     = ptrace.(solve(problemA, Nouts = [6], reltol = 1e-10)[2], 1)

  problemB = WaveguideProblem((σm+σp, [], σm, spindown(testBasis)),
		WavePacket(SoftBoxMode(τ = 1.0), Coherent(1.0)),
		SoftBoxMode(τ = 2.0), 3.
	)
  solB     = solve(problemB, Nouts = [6], reltol = 1e-10)[2]

  @test sum(tracedistance.(solA, solB)) < 1e-4
end

@testset "ContinuesWave(::Function) - normalisable" begin
	# ----------------------------------------
	# Show that ContinuesWave(::Function) creates correct state
	testBasis = SpinBasis(1//2)
	σm, σp = sigmam(testBasis), sigmap(testBasis)

	ts = range(-5.0, 5.0, length = 101)
	problemA = WaveguideProblem((σm+σp, [], σm, spindown(testBasis)),
    WavePacket(GaussMode(), Coherent( (2π)^(1/4) )), ts
  )
	solA = solve(problemA, Nouts = [6], reltol = 1e-10)[2]

	problemB = WaveguideProblem((σm+σp, [], σm, spindown(testBasis)),
    ContinuousWave(t -> e^(-t^2 / 2)), ts
  )
	solB = solve(problemA, Nouts = [6], reltol = 1e-10)[2]

	@test sum(tracedistance.(solA, solB)) ≈ 0.0
end

@testset "Controlling the output Fock space" begin
	# ----------------------------------------
	# Setting the cutoff vs. setting `Nouts`
  testBasis = SpinBasis(1//2)
  σm, σp = sigmam(testBasis), sigmap(testBasis)

  problemA = WaveguideProblem((σm+σp, [], σm, spindown(testBasis)),
    WavePacket(SoftBoxMode(τ = 1.0), Coherent(1.0, 10)),
		SoftBoxMode(τ = 2.0), 3.
  )
  solA     = solve(problemA, reltol = 1e-10)[2]

  problemB = WaveguideProblem((σm+σp, [], σm, spindown(testBasis)),
		WavePacket(SoftBoxMode(τ = 1.0), Coherent(1.0)),
		SoftBoxMode(τ = 2.0), 3.
	)
  solB     = solve(problemB, Nouts = [10], reltol = 1e-10)[2]

  @test sum(tracedistance.(solA, solB)) < 1e-10
end

@testset "DrivenProblem - DisplacedFock" begin
	# ----------------------------------------
	# Compare the `DisplacedFock` implementation to the manually created version
  testBasis = SpinBasis(1//2)
  σm, σp = sigmam(testBasis), sigmap(testBasis)

	aᵀ, ψα = create(FockBasis(13)), coherentstate(FockBasis(13), 1.0)
	data = (aᵀ * ψα - ψα).data
  problemA = WaveguideProblem((σm+σp, [], σm, spindown(testBasis)),
    WavePacket(SoftBoxMode(τ = 1.0), ArbitraryState(data)), 3.
  )
  solA     = ptrace.(solve(problemA, reltol = 1e-10)[2], 1)

  problemB = WaveguideProblem((σm+σp, [], σm, spindown(testBasis)),
		WavePacket(SoftBoxMode(τ = 1.0), DisplacedFock(1.0, 1)), 3.
	)
  solB     = ptrace.(solve(problemB, reltol = 1e-10)[2], 1)

  @test sum(tracedistance.(solA, solB)) < 1e-5
end

@testset "DrivenProblem - DisplacedArbitraryState" begin
	# ----------------------------------------
	# Compare the `DisplacedArbitraryState` implementation to the manually created version
  testBasis = SpinBasis(1//2)
  σm, σp = sigmam(testBasis), sigmap(testBasis)

	aᵀ, ψα = create(FockBasis(13)), coherentstate(FockBasis(13), 1.0)
	data = (aᵀ * ψα - ψα).data
  problemA = WaveguideProblem((σm+σp, [], σm, spindown(testBasis)),
    WavePacket(SoftBoxMode(τ = 1.0), ArbitraryState(data)), 3.
  )
  solA     = ptrace.(solve(problemA, reltol = 1e-10)[2], 1)

  problemB = WaveguideProblem((σm+σp, [], σm, spindown(testBasis)),
		WavePacket(SoftBoxMode(τ = 1.0), DisplacedArbitraryState(1.0, [0.0im, 1.0])), 3.
	)
  solB     = ptrace.(solve(problemB, reltol = 1e-10)[2], 1)

  @test sum(tracedistance.(solA, solB)) < 1e-5
end

@testset "ScatterProblem - DisplacedFock" begin
	# ----------------------------------------
	# Compare the `DisplacedFock` implementation to the manually created version
  testBasis = SpinBasis(1//2)
  σm, σp = sigmam(testBasis), sigmap(testBasis)

	aᵀ, ψα = create(FockBasis(13)), coherentstate(FockBasis(13), 1.0)
	data = (aᵀ * ψα - ψα).data
  problemA = WaveguideProblem((σm+σp, [], σm, spindown(testBasis)),
    WavePacket(SoftBoxMode(τ = 1.0), ArbitraryState(data)), SoftBoxMode(τ = 2.0), 3.
  )
  solA     = ptrace.(solve(problemA, reltol = 1e-10)[2], 1)

  problemB = WaveguideProblem((σm+σp, [], σm, spindown(testBasis)),
		WavePacket(SoftBoxMode(τ = 1.0), DisplacedFock(1.0, 1)), SoftBoxMode(τ = 2.0), 3.
	)
  solB     = ptrace.(solve(problemB, Nouts = [13], reltol = 1e-10)[2], 1)

  @test sum(tracedistance.(solA, solB)) < 1e-4
end

@testset "DrivenProblem - DisplacedArbitraryState" begin
	# ----------------------------------------
	# Compare the `DisplacedArbitraryState` implementation to the manually created version
  testBasis = SpinBasis(1//2)
  σm, σp = sigmam(testBasis), sigmap(testBasis)

	aᵀ, ψα = create(FockBasis(13)), coherentstate(FockBasis(13), 1.0)
	data = (aᵀ * ψα - ψα).data
  problemA = WaveguideProblem((σm+σp, [], σm, spindown(testBasis)),
    WavePacket(SoftBoxMode(τ = 1.0), ArbitraryState(data)), SoftBoxMode(τ = 2.0), 3.
  )
  solA     = ptrace.(solve(problemA, reltol = 1e-10)[2], 1)

  problemB = WaveguideProblem((σm+σp, [], σm, spindown(testBasis)),
		WavePacket(SoftBoxMode(τ = 1.0), DisplacedArbitraryState(1.0, [0.0im, 1.0])),
		SoftBoxMode(τ = 2.0), 3.
	)
  solB     = ptrace.(solve(problemB, Nouts = [13], reltol = 1e-10)[2], 1)

  @test sum(tracedistance.(solA, solB)) < 1e-4
end
