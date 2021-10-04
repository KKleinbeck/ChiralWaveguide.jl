include("./Compressions.jl")
export Algebraic, Exponential, Trigonometric

include("./Couplings.jl")
export generateCoupling, modeContribution, modifiedMode

include("./ModeFunctions.jl")
export Mode, HardBoxMode, SoftBoxMode, SoftBoxExpMode, GaussMode, ExpMode

include("./WavePacketStates.jl")
export WavePacketState, Displaced, NonDisplaced,
  Coherent, ArbitraryState, Fock, SqueezedVacuum, createState,
  DisplacedArbitraryState, DisplacedFock

abstract type _WavePacket{T <: WavePacketState} end

mutable struct WavePacket{S} <: _WavePacket{S}
  mode::Mode
  state::S
end

mutable struct ContinuousWave <: _WavePacket{Coherent}
  α::Number
end

function ContinuousWave(α::Function)
	warning = "`ContinuousWave(α::Function)` doesn't define gᵢ, gₒ, or norm, " *
		"use `WavePacket(Mode(α), Coherent(1.0))` instead"
  mode = Mode(α, # Since we do not want to assume `α` to be normalised.
    t -> begin @warn(warning); NaN end,
    t -> begin @warn(warning); NaN end,
		t -> begin @warn(warning); NaN end
  )
  WavePacket(mode, Coherent(1.0))
end

NoPacket() = ContinuousWave(0.0)
