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
		FourierMode(1/√(5), [1/√(5), 1/√(5)], [1/√(5), 1/√(5)], σ = 3.0),
		HardBoxMode(σ = 3.0),    # give a width to help the integrator
		SoftBoxMode(σ = 3.0),    # A little bit of fine tuning is involved here, but
		SoftBoxExpMode(σ = 2.5), # I simple want to test the normalisation conditions.
		GaussMode(),
		ExpMode(γ =  4.0),
		ExpMode(γ = -4.0)
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

@testset "FourierMode and HardBoxMode equivalence" begin
	hb = HardBoxMode()
	f  = FourierMode(1., [], [])

	ts = range(0.0, 1.0, length = 101)

	@test sum(abs, hb.modeFunction.(ts) .- f.modeFunction.(ts)) / length(ts) < 1e-9
	@test sum(abs, hb.gᵢ.(ts) .- f.gᵢ.(ts)) / length(ts) < 1e-9
	@test sum(abs, hb.gₒ.(ts) .- f.gₒ.(ts)) / length(ts) < 1e-9
	@test sum(abs, hb.norm.(ts) .- f.norm.(ts)) / length(ts) < 1e-9
end

using QuantumOptics
@testset "Cavity Hopping" begin
	# Idea of this test is, that a cavity for a given mode
	# should have 100% into a cavity for the same mode.
	modes = [
		FourierMode(rand(), rand(2), rand(2), t₀ = -5.0),
		HardBoxMode(t₀ = -5.0),
		SoftBoxMode(),
		SoftBoxExpMode(),
		GaussMode(),
		ExpMode(γ =  3.5), # Needed to fit in the time window
		ExpMode(γ = -3.5)
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
		  ts, ρs = timeevolution.master_dynamic(
				[-5:0.01:5.0;], ψ₀, (t,ψ) -> (H(t), [L(t)], [L(t)']), reltol = 1e-7, abstol = 1e-9
			)

		  @test abs(expect(2, a'a, ρs[end]) - n) < 1e-7
		end
	end
end
