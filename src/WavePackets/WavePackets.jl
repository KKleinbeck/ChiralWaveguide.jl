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

NoPacket() = ContinuousWave(0.0)
