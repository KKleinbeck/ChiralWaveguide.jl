"""
    WavePacketState

Abstract base class for all states.
"""
abstract type WavePacketState end

"""
    Coherent(α[, N_cutoff])

Describes a coherent state ``|α⟩``. The autmatic choice for `N_cutoff`
yields a normalisation of above 0.999 for α ∈ [0,10].
"""
struct Coherent <: WavePacketState
	α::Complex{Float64}
	N_cutoff::Int
end

coherent_cutoff(α)::Int = ceil(Int, abs(α) * (abs(α) + 3.5) + 0.45 * √(abs(α)))
Coherent(α)                    = Coherent(Complex{Float64}(α), coherent_cutoff(α))
Coherent(α::Float64, N_cutoff) = Coherent(Complex{Float64}(α), N_cutoff)

abstract type Displaced    <: WavePacketState end
abstract type NonDisplaced <: WavePacketState end

#-------------------------------------------------------
# NonDisplaced States

"""
    ArbitraryState(amplitudes[, N_cutoff])

Creates a state with specific amplitudes.
"""
struct ArbitraryState <: NonDisplaced
	amplitudes::Array{Complex{Float64}, 1}
	N_cutoff::Int
end
ArbitraryState(amps::Array{Complex{Float64}, 1}) = ArbitraryState(amps, length(amps) - 1)

"""
    Fock(n[, N_cutoff])

Decribes the Fock state ``|n⟩``.
"""
struct Fock <: NonDisplaced
	n::Int
	N_cutoff::Int
end
Fock(n) = Fock(n, n)

"""
    SqueezedVacuum(r, ϕ[, N_cutoff])

Decribes the a squeezed state with squeezing amplitude `r` and squeezing angle `ϕ`, i.e., the
state ``\\exp[(ξ^{*} a^2 - ξ a^{\\dagger 2})/2] |0⟩`` with ``ξ = r e^{i ϕ}``.
"""
struct SqueezedVacuum <: NonDisplaced
	r::Float64
	ϕ::Float64
	N_cutoff::Int
end
# magic number choosen so that for r ∈ [0,3] the normalisation of the state is above 0.99
SqueezedVacuum(r, ϕ = 0.0) = SqueezedVacuum(r, ϕ, 2 * ceil(Int, -1.7 / log(tanh(r))))
SqueezedVacuum(ξ::Complex{Float64}) = SqueezedVacuum(abs(ξ), angle(ξ))
SqueezedVacuum(ξ::Complex{Float64}, N_cutoff) = SqueezedVacuum(abs(ξ), angle(ξ), N_cutoff)

#-------------------------------------------------------
# Displaced States

"""
    DisplacedArbitraryState(α, amplitudes[, N_cutoff])

Describes a state with specific amplitudes, which then is displaced by ``D(α)``.
"""
struct DisplacedArbitraryState <: Displaced
	α::Complex{Float64}
	amplitudes::Array{Complex{Float64}, 1}
	N_cutoff::Int
end
DisplacedArbitraryState(α, amps::Array{Complex{Float64}, 1}) =
	DisplacedArbitraryState(Complex{Float64}(α), amps, length(amps) - 1)

"""
	  Fock(n[, N_cutoff])

Decribes the displaced Fock state ``D(α)|n⟩``.
"""
struct DisplacedFock <: Displaced
	α::Complex{Float64}
	n::Int
	N_cutoff::Int
end
DisplacedFock(α, n) = DisplacedFock(Complex{Float64}(α), n, n)

#-------------------------------------------------------
# Implementation of `createState`

createState(wpt::Union{ArbitraryState, DisplacedArbitraryState}) = wpt.amplitudes

function createState(wpt::Union{Fock, DisplacedFock})
	data = zeros(wpt.N_cutoff + 1)
	wpt.n ≤ wpt.N_cutoff && (data[wpt.n+1] = 1.0)
	return data
end

function createState(wpt::SqueezedVacuum)
	data = zeros(Complex{Float64}, wpt.N_cutoff + 1)

	r, ϕ = wpt.r, wpt.ϕ
	tanh_r, exp_i_ϕ = tanh(r), exp(1im * ϕ)
	for n ∈ 0:wpt.N_cutoff÷2
		data[2n+1] = ( -exp_i_ϕ * tanh_r )^n *
			( √(factorial(2.0*n)) / (2^n * factorial(1.0*n)) )
	end
	return data ./ √(sum(abs2, data))
end
