using DifferentialEquations: ODEProblem, DifferentialEquations

@inline function __integrate(u::Function, domain::Tuple{T,T}; kwargs...) where T <: AbstractFloat
	prob = ODEProblem( (f, p, τ) -> u(τ), zero(T), domain)
	return DifferentialEquations.solve(prob; kwargs...)
end

@inline function __getNorm(u::Function, c::Compression; kwargs...)
	σ, μ = c.σ, c.μ
	sol = __integrate(
		τ -> abs(τ) ≈ 1.0 ? 0.0 : σ * abs2( u(σ * c.decomp(τ) + μ) ) * c.inv_jacobian(τ),
		(-1.0, 1.0);
		kwargs...
	)
	return t -> abs(sol(c.comp((t - μ) / σ))), sol(1.0)
end

"""
    generateCouplings(u; compression, kwargs...)

inputs
- u:            mode function
- compression:  which compression of ℝ → (-1, 1) will be used for the differential equation
								Options are the keys of the `Compressions` dictionary
								or any `Compression` Object
- kwargs:				all other kwargs get's passed to the numerical integrator

output:
- (gᵢ, gₒ, norm): the coupling constant calculated from the mode function and its norm function.
"""
# TODO test with (df, f, p, τ) lambda expression

function generateCoupling(u::Function; compression::Union{Symbol, Compression} = :exponential, kwargs...)
	if isa(compression, Symbol)
		@assert compression ∈ keys(Compressions)

		return _generateCoupling(u, Compressions[compression]; kwargs...)
	end

	return _generateCoupling(u, compression; kwargs...)
end

function _generateCoupling(u::Function, c::Compression; kwargs...)
	norm, Norm = __getNorm(u, c; kwargs...)
	gᵢ(t) =   u(t) / √(abs(Norm - norm(t)) + (u(t) == 0.0))
	gₒ(t) = - u(t) / √(norm(t) + (u(t) == 0.0))

	return gᵢ, gₒ, norm
end

# """
#     modeContribution(u, v; compression)
#
# inputs
# - u:            mode function of the cavity
# - v:						mode function falsely contributing to the cavity
# - compression:  which compression of ℝ → (-1, 1) will be used for the differential equation
# 								Options are the keys of the `Compressions` dictionary
# 								or any `Compression` Object
#
# output:
# - α:   the contribution of mode u to the cavity
# """
# # TODO test with (df, f, p, τ) lambda expression
# function modeContribution(u::Mode, v::Mode; compression::Union{Symbol, Compression} = :exponential)
# 	if isa(compression, Symbol)
# 		@assert compression ∈ keys(Compressions)
#
# 		return _modeContribution(u, v, Compressions[compression])
# 	end
#
# 	return _modeContribution(u, v, compression)
# end
#
# function _modeContribution(uM::Mode, vM::Mode, c::Compression)
# 	u, v = uM.modeFunction, vM.modeFunction
# 	dfdτ(f, p, τ) = begin
# 		abs(τ) == 1.0 && return 0.0
#
# 		t = c.decomp(τ)
# 		df = u(t) * v(t) * c.inv_jacobian(τ)
# 	end
#
# 	prob = ODEProblem(dfdτ, 0., (-1.0, 1.0))
#   sol = solve(prob, reltol = 1e-5, abstol = 1e-8)
#
# 	return t -> sol(c.comp(t)) / √(uM.norm(t))
# end
#
# """
#     modifiedMode(u, v; compression)
#
# inputs
# - u:            `Mode` of the first cavity
# - v:						target `Mode` for the second cavity
# - [α]:					predetermined modeContribution of v to the first cavity. Optional.
# 								Useful, if α should be calculated with a different compression than the
# 								modified mode.
# - compression:  which compression of ℝ → (-1, 1) will be used for the differential equation
# 								Options are the keys of the `Compressions` dictionary
# 								or any `Compression` Object
#
# output:
# - ṽ:   modified `Mode` for the second cavity, accounting for the scattering at the first cavity
# """
# # TODO test with (df, f, p, τ) lambda expression
# function modifiedMode(u::Mode, v::Mode, α::Union{Function, Nothing} = nothing;
# 		compression::Union{Symbol, Compression} = :exponential
# 	)
#
# 	if isa(compression, Symbol)
# 		@assert compression ∈ keys(Compressions)
# 		compression = Compressions[compression]
# 	end
#
# 	isnothing(α) && (α = _modeContribution(u, v, compression))
# 	return _modifiedMode(u, v, α, compression)
# end
#
# function _modifiedMode(u::Mode, v::Mode, α::Function, compression::Compression)
# 	vMod(t) = v.modeFunction(t) + u.gₒ(t) * α(t)
#
# 	return Mode(vMod, compression = compression)
# end
