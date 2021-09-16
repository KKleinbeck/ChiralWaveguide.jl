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

  @test tracedistance.(solA, solB) |> sum < 1e-5
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

  @test tracedistance.(solA, solB) |> sum < 1e-4
end
