abstract type WavePacketState end

struct Coherent <: WavePacketState
	α::Complex{Float64}
	N_cutoff::Int
end
# magic numbers choosen so that for α ∈ [0,10] the normalisation of the state is above 0.999
coherent_cutoff(α)::Int = ceil(Int, abs(α) * (abs(α) + 3.5) + 0.45 * √(abs(α)))
Coherent(α)                    = Coherent(Complex{Float64}(α), coherent_cutoff(α))
Coherent(α::Float64, N_cutoff) = Coherent(Complex{Float64}(α), N_cutoff)

abstract type Displaced    <: WavePacketState end
abstract type NonDisplaced <: WavePacketState end

#-------------------------------------------------------
# NonDisplaced States

struct ArbitraryState <: NonDisplaced
	data::Array{Complex{Float64}, 1}
	N_cutoff::Int
end
ArbitraryState(data::Array{Complex{Float64}, 1}) = ArbitraryState(data, length(data) - 1)
createState(wpt::ArbitraryState) = wpt.data

struct Fock <: NonDisplaced
	n::Int
	N_cutoff::Int
end
Fock(n) = Fock(n, n)
function createState(wpt::Fock)
	data = zeros(wpt.N_cutoff + 1)
	if wpt.n ≤ wpt.N_cutoff
		data[wpt.n+1] = 1.0
	end
	return data
end

struct SqueezedVacuum <: NonDisplaced
	r::Float64
	ϕ::Float64
	N_cutoff::Int
end
# magic number choosen so that for r ∈ [0,3] the normalisation of the state is above 0.99
SqueezedVacuum(r, ϕ = 0.0) = SqueezedVacuum(r, ϕ, 2 * ceil(Int, -1.7 / log(tanh(r))))
SqueezedVacuum(ξ::Complex{Float64}) = SqueezedVacuum(abs(ξ), angle(ξ))
SqueezedVacuum(ξ::Complex{Float64}, N_cutoff) = SqueezedVacuum(abs(ξ), angle(ξ), N_cutoff)
function createState(wpt::SqueezedVacuum)
	data = zeros(wpt.N_cutoff + 1)

	r, ϕ = wpt.r, wpt.ϕ
	tanh_r, exp_i_ϕ = tanh(r), exp(1im * ϕ)
	for n ∈ 0:wpt.N_cutoff÷2
		data[2n+1] = ( -exp_i_ϕ * tanh_r )^n *
			(√(factorial(2.0*n)) / (2^n * factorial(1.0*n)) )
	end
	return data ./ √(sum(abs2, data))
end

#-------------------------------------------------------
# Displaced States

struct DisplacedArbitraryState <: Displaced
	α::Complex{Float64}
	data::Array{Complex{Float64}, 1}
	N_cutoff::Int
end
DisplacedArbitraryState(α, data::Array{Complex{Float64}, 1}) =
	DisplacedArbitraryState(Complex{Float64}(α), data, length(data) - 1)
createState(wpt::DisplacedArbitraryState) = wpt.data

struct DisplacedFock <: Displaced
	α::Complex{Float64}
	n::Int
	N_cutoff::Int
end
DisplacedFock(α, n) = DisplacedFock(Complex{Float64}(α), n, n)
function createState(wpt::DisplacedFock)
	data = zeros(wpt.N_cutoff + 1)
	if wpt.n ≤ wpt.N_cutoff
		data[wpt.n+1] = 1.0
	end
	return data
end
