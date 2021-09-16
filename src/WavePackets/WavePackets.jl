include("./Compressions.jl")

include("./Couplings.jl")
export generateCoupling, modeContribution, modifiedMode

include("./ModeFunctions.jl")
export Mode, HardBoxMode, SoftBoxMode, SoftBoxExpMode, GaussMode, FlatMode

include("./WavePacketStates.jl")
export WavePacketState, Displaced, NonDisplaced,
  Coherent, ArbitraryState, Fock, SqueezedVacuum, createState

mutable struct WavePacket{T <: WavePacketState}
  mode::Mode
  state::T
end
