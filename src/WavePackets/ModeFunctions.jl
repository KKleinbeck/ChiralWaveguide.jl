using SpecialFunctions, HypergeometricFunctions

"""
    Mode(modeFunction; compression = :algebraic, kwargs...)
    Mode(modeFunction, gᵢ, gₒ, norm)

Container for a cavity mode.

# Arguments
- `modeFunction`: the wave function of the mode
- `gᵢ`: coupling function for an input cavity
- `gₒ`: coupling function for an output cavity
- `norm`: square integral of the mode up to time t

If only the mode function is provided, then all other quantities are numerically determined.
This is done by mapping ``ℝ → (-1, 1)`` and then numerically solving the resulting integrals.
The argument `compression` picks the map ``ℝ → (-1, 1)`` and can be a [`Compression`](@ref) or
the respective symbol.
Allowed symbols are the keys of the `ChiralWaveguide.Compressions` dictionary.
Additional keyword arguments will be passed to the numerical integrator.
"""
mutable struct Mode
	modeFunction::Function
	gᵢ::Function # Input coupling
	gₒ::Function # Output coupling
	norm::Function # Square integral of the mode function
end

function Mode(mF::Function; compression::Union{Symbol, Compression} = :algebraic, kwargs...)
	gᵢ, gₒ, norm = generateCoupling(mF; compression = compression, kwargs...)
	Mode(mF, gᵢ, gₒ, norm)
end

#-------------------------------------------------------
# Special predefined and numerically accurate modes


"""
    HardBoxMode(; t₀ = 0.0, σ = 1.0)

The mode is a box function, i.e., the mode takes the value ``1 / \\sqrt{σ}`` for ``t₀ < t < t₀ + σ``
and 0 otherwise.
"""
function HardBoxMode(; t₀ = 0.0, σ = 1.0)
	function norm(t)
		t < t₀     && return 0.0
		t < t₀ + σ && return t / σ
		return 1.0
	end
  return Mode(
    t -> (t₀ < t) && (t < t₀ + σ) ?  1.0 / √(σ)          : 0.0,
    t -> (t₀ < t) && (t < t₀ + σ) ?  1.0 / √(σ - t + t₀) : 0.0,
    t -> (t₀ < t) && (t < t₀ + σ) ? -1.0 / √(    t - t₀) : 0.0,
		norm
  )
end



"""
    SoftBoxMode(; τ = 0.0, σ = 1.0, n::Int = 10)

The mode is an approximation to the box mode,
``\\mathrm{box}(t) ≈ \\frac{1}{\\sqrt{(2t)^{2n} + 1}}``.
The parameters `τ` and `σ` control the center and width respectively, `n` the exponent.
"""
function SoftBoxMode(; τ = 0.0, σ = 1.0, n::Int = 10)
  n̄ = 2n
  norm = π * csc(π / n̄) / n̄

	function softBox(t; n::Int = 20)
	  norm = π * csc(π / 2n) / 2n
	  return (1 / √(norm) ) / √( (2t)^(2n) + 1 )
	end

	function softBoxNorm(t)
		return t * _₂F₁(1, 1/n̄, 1.0 + 1/n̄, -(2t)^n̄) / norm + 1/2
	end

  function softBoxCoupling(t)
    if abs(t) < 0.5
      return softBox(t, n = n) / √(1.0 - softBoxNorm(t))
    end
    Θmt = t < 0.0 ? 1.0 : 0.0
    twotⁿ = (2t)^(n̄)

    mainTerm = (-t / (1 - n̄)) * _₂F₁(1, 1, 2 - 1/n̄, 1/(1.0 + twotⁿ) )

    return 1.0 / √(
      mainTerm + norm * (1.0 + twotⁿ) * Θmt
    )
	end

  return Mode(
    t ->  softBox( (t - τ) / σ, n = n)  / √(σ),
    t ->  softBoxCoupling( (t - τ) / σ) / √(σ),
    t -> -softBoxCoupling(-(t - τ) / σ) / √(σ),
		t ->  softBoxNorm(     (t - τ) / σ) / √(σ)
  )
end


"""
    SoftBoxExpMode(; τ = 0.0, σ = 1.0, γ = 10.0)

A box mode with exponential decaying flanks, given by

                       ⌜ 1                     for |t| < 1/2
    softBoxExp(t, γ) = |
                       ⌞ exp(-γ(|t| - 1/2))    for |t| > 1/2

The parameters `τ` and `σ` control the center and width respectively, `γ` the deay rate.
"""
function SoftBoxExpMode(; τ = 0.0, σ = 1.0, γ = 10.0)
	function softBoxExp(t; γ = 10.0)
	  t = abs(t)
	  return t > 0.5 ? exp(-γ * (t-1/2)) : 1.0
	end

  function softBoxExpNorm(t)
    if t < -0.5
      return exp(2γ * (t + 1/2)) / 2γ
    elseif t < 0.5
      return 1/(2γ) + t + 1/2
    end

    return 1/γ + 1.0 - exp(-2γ * (t - 1/2)) / 2γ
	end

  function softBoxExpCoupling(t)
    if t < -0.5
      return exp(γ * (t+1/2)) / √(1 + (1/γ) - exp(2γ * (t+1/2))/(2γ))
    elseif t < 0.5
      return 1/√(1/2 + 1/(2γ) - t)
    end

    return √(2γ)
  end

  return Mode(
    t ->  √(γ / (γ * σ + σ)) * softBoxExp( (t - τ) / σ, γ = γ),
    t ->  softBoxExpCoupling( (t - τ) / σ) / √(σ),
    t -> -softBoxExpCoupling(-(t - τ) / σ) / √(σ),
		t ->  √(γ / (γ * σ + σ)) * softBoxExpNorm( (t - τ) / σ)
  )
end


"""
    GaussMode(; τ = 0.0, σ = 1.0)

Numerical stable mode function & couplings of Gaussian wave packet with mean τ and variance σ.
"""
function GaussMode(; τ = 0.0, σ = 1.0)
  return Mode(
    t ->  exp(-(t - τ)^2 / (2 * σ^2) ) / √( σ * √(π)),
    t ->  exp(-(t - τ)^2 / (2 * σ^2) ) / √((σ * √(π)) * erfc( (t - τ) / σ) / 2 + eps(0.0)),
    t -> -exp(-(t - τ)^2 / (2 * σ^2) ) / √((σ * √(π)) * erfc(-(t - τ) / σ) / 2 + eps(0.0)),
		t ->  (1 + erf( (t - τ) / σ)) / 2.0
  )
end


"""
    ExpMode(; t₀ = 0.0, γ = 1.0)

Exponentially raises up to time t₀ or decays after t₀, depending on the sign of the rate γ.
Amplitude decays with rate γ/2, so that the probability density decays with γ.
"""
function ExpMode(; t₀ = 0.0, γ = 1.0)
	@assert γ != 0.0

	if γ < 0.0
		return Mode(
			t -> t > t₀ ?  √(-γ) * exp(γ*(t-t₀)/2) : 0,
			t -> t > t₀ ?  √(-γ) : 0.0,
			t -> t > t₀ ? -√(-γ) * exp(γ*(t-t₀)/2) / √( 1.0 - exp(γ*(t-t₀)) ) : 0.0,
			t -> t > t₀ ? 1.0 - exp(γ*(t-t₀)) : 0.0
		)
	end
	Mode(
		t -> t < t₀ ?  √(γ) * exp(γ*(t-t₀)/2) : 0,
		t -> t < t₀ ?  √(γ) * exp(γ*(t-t₀)/2) / √( 1.0 - exp(γ*(t-t₀)) ) : 0.0,
		t -> t < t₀ ? -√(γ) : 0.0,
		t -> t < t₀ ? exp(γ*(t-t₀)) : 1.0
	)
end
