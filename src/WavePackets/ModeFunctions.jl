using SpecialFunctions, HypergeometricFunctions

"""
Mode

Container for a cavity mode. Consists of:
- `modeFunction`: the wave function of the mode
- `gᵢ`: coupling function for an input cavity
- `gₒ`: coupling function for an output cavity
- `norm`: square integral of the mode up to time t

If only the mode function is provided, then all other
quantities are numerically determined. This is done by
first mapping ℝ → (-1, 1) and then numerically solving
the resulting integrals. For the map ℝ → (-1, 1) the
compressions may be picked by using
`Mode(modeFunction, compression::Union{Symbol, Compression})`.
Allowed Symbols are the keys of the `Compressions` dictionary.
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
HardBoxMode(; t₀, σ)

Gives the mode unit box function, i.e., the mode takes values
	1 / √(σ)   for t₀ < t < t₀ + σ
	0          otherwise
"""
function HardBoxMode(; t₀ = 0.0, σ = 1.0)
	function norm(t)
		x < t₀     && return 0.0
		x < t₀ + σ && return t / σ
		return 1.0
	end
  return Mode(
    x -> (t₀ < x) && (x < t₀ + σ) ?  1.0 / √(σ)          : 0.0,
    x -> (t₀ < x) && (x < t₀ + σ) ?  1.0 / √(σ - x + t₀) : 0.0,
    x -> (t₀ < x) && (x < t₀ + σ) ? -1.0 / √(    x - t₀) : 0.0,
		norm
  )
end


"""
softBox(x, n)

The UnitBox may be approximated by 1/(x^2n + 1). However, this function is hard
for exact calculations and a better alternative is (which is returned normalized):
	1 / sqrt(x^2n + 1)
"""
function softBox(x::T; n::Int = 20) where {T<:AbstractFloat}
  norm = π * csc(π / 2n) / 2n
  return (one(T) / √(norm) ) / √( (2x)^(2n) + one(T) )
end

function SoftBoxMode(; τ = 0.0, σ = 1.0, n::Int = 10)
  n̄ = 2n
  norm = π * csc(π / n̄) / n̄

	function softBoxNorm(x)
		return x * _₂F₁(1, 1/n̄, 1.0 + 1/n̄, -(2x)^n̄) / norm + 1/2
	end

  function softBoxCoupling(x)
    if abs(x) < 0.5
      return softBox(x, n = n) / √(1.0 - softBoxNorm(x))
    end
    Θmx = x < 0.0 ? 1.0 : 0.0
    twoxⁿ = (2x)^(n̄)

    mainTerm = (-x / (1 - n̄)) * _₂F₁(1, 1, 2 - 1/n̄, 1/(1.0 + twoxⁿ) )

    return 1.0 / √(
      mainTerm + norm * (1.0 + twoxⁿ) * Θmx
    )
	end

  return Mode(
    x ->  softBox( (x - τ) / σ, n = n)  / √(σ),
    x ->  softBoxCoupling( (x - τ) / σ) / √(σ),
    x -> -softBoxCoupling(-(x - τ) / σ) / √(σ),
		x ->  softBoxNorm( (x - τ) / σ) / √(σ)
  )
end


"""
softBoxExp(x, γ)

Approximate the UnitBox by a fast decaying and continous function. given by

                       ⌜ 1                     for |x| < 1/2
    softBoxExp(x, γ) = |
                       ⌞ exp(-γ(|x| - 1/2))    for |x| > 1/2
"""
function softBoxExp(x; γ = 10.0)
  x = abs(x)
  return x > 0.5 ? exp(-γ * (x-1/2)) : one(x)
end

function SoftBoxExpMode(; τ = 0.0, σ = 1.0, γ = 10.0)
  function softBoxExpNorm(x)
    if x < -0.5
      return exp(2γ * (x + 1/2)) / 2γ
    elseif x < 0.5
      return 1/(2γ) + x + 1/2
    end

    return 1/γ + 1.0 - exp(-2γ * (x - 1/2)) / 2γ
	end

  function softBoxExpCoupling(x)
    if x < -0.5
      return exp(γ * (x+1/2)) / √(1 + (1/γ) - exp(2γ * (x+1/2))/(2γ))
    elseif x < 0.5
      return 1/√(1/2 + 1/(2γ) - x)
    end

    return √(2γ)
  end

  return Mode(
    x ->  √(γ / (γ * σ + σ)) * softBoxExp( (x - τ) / σ, γ = γ),
    x ->  softBoxExpCoupling( (x - τ) / σ) / √(σ),
    x -> -softBoxExpCoupling(-(x - τ) / σ) / √(σ),
		x ->  √(γ / (γ * σ + σ)) * softBoxExpNorm( (x - τ) / σ)
  )
end


"""
  GaussMode(τ, σ)

  Numerical stable mode function & couplings of Gaussian wave packet with mean τ and variance σ.
"""
function GaussMode(; τ = 0.0, σ = 1.0)
  return Mode(
    x ->  exp(-(x - τ)^2 / (2 * σ^2) ) / √( σ * √(π)),
    x ->  exp(-(x - τ)^2 / (2 * σ^2) ) / √((σ * √(π)) * erfc( (x - τ) / σ) / 2 + eps(0.0)),
    x -> -exp(-(x - τ)^2 / (2 * σ^2) ) / √((σ * √(π)) * erfc(-(x - τ) / σ) / 2 + eps(0.0)),
		x ->  (1 + erf( (x - τ) / σ)) / 2.0
  )
end


"""
FlatMode(tf)

Constant wave, without any couplings. `tf` gives the final time, i.e., a cutoff.
"""
function FlatMode(; tf = Inf)
  return Mode(
    x -> x < tf ? 1.0 : 0.0,
    x -> begin @warn("`flatMode` doesn't define gₒ"); NaN end,
    x -> begin @warn("`flatMode` doesn't define gᵢ"); NaN end,
		x -> begin @warn("`flatMode` doesn't define a norm"); NaN end
  )
end
