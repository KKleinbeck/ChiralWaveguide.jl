using QuadGK

@testset "Integration and Coupling generation" begin
	# ----------------------------------------
	# Test the numerical integration capabilities
	antiDerX = ChiralWaveguide.__integrate(x -> x, (0.0, 1.0))

	@test antiDerX(1.0) ≈ 1/2	# local accuracy
	@test sum(x -> abs(antiDerX(x) - x^2 / 2), [0.0:0.01:1.0;]) ≤ 1e-12 # global accuracy
end

@testset "Modes with offset" begin
	μ, σ = 10.0, 0.1

	# ----------------------------------------
	# Test normalisation
	compressions = [Algebraic(σ = σ, μ = μ), Exponential(σ = σ, μ = μ), Trigonometric(σ = σ, μ = μ)]
	for compression ∈ compressions
		testMode = Mode(
			t -> abs(t - μ) < σ/2 ? √(1/σ) : 0.0, compression = compression,
			abstol = 1e-8, reltol = 1e-10
		)
		@test abs(testMode.norm(2μ) - 1.0) < 1e-7
	end
end

@testset "Numerical Mode Couplings" begin
	modes = [
		GaussMode(),
		HardBoxMode(σ = 3.0), # since quadgk otherwise doesn't hit a non zero point
		SoftBoxMode(σ = 3.0), # due to the slow decay
		SoftBoxExpMode(σ = 2.5)
	]
	# ----------------------------------------
	# Test Numerical Derivations
	ts = range(-3.0, 3.0, length = 100)
	for mode ∈ modes
		numericalMode = Mode(t -> mode.modeFunction(t), reltol = 1e-8, abstol = 1e-12)
		@test sum(abs, numericalMode.gₒ.(ts) .- mode.gₒ.(ts)) / length(ts) < 1e-7
	end

	# ----------------------------------------
	# Integration test
	for mode ∈ modes
		@test quadgk(
			t -> mode.modeFunction(t) * mode.gₒ(t),
			-10.0, 10.0
		)[1] ≈ -2.0
	end
end

using QuantumOptics
@testset "Cavity Hopping" begin
	# Idea of this test is, that a cavity for a given mode
	# should have 100% into a cavity for the same mode.
	modes = [
		GaussMode(),
		HardBoxMode(t₀ = -5.0),
		SoftBoxMode(),
		SoftBoxExpMode()
	]

	for mode ∈ modes
		for n ∈ [1, 2, 3] # Fock state scattering
		  basis = FockBasis(n)
		  a = destroy(basis)

		  gᵢ, gₒ = mode.gᵢ, mode.gₒ

		  H(t) = 0.5im * gᵢ(t) * gₒ(t) * (a' ⊗ a - a ⊗ a')
		  L(t) = gᵢ(t) * a ⊗ one(basis) + gₒ(t) * one(basis) ⊗ a

		  ψ₀ = fockstate(basis, n) ⊗ fockstate(basis, 0)

		  local ts, ρs
		  ts, ρs = timeevolution.master_dynamic([-5:0.01:5.0;], ψ₀, (t,ψ) -> (H(t), [L(t)], [L(t)']))

		  @test abs(expect(2, a' * a, ρs[end]) - n) < 1e-7
		end
	end
end
