mutable struct _ScatterProblem{O, S <: WavePacketState} <: WaveguideProblem
	H::O
	Ls::Array{O}
	σ::O
	system_state

	ψᵢ::WavePacket{S}
	ψₒ::Array{Mode}

	ts::StepRangeLen
end

function WaveguideProblem(sys, ψᵢ::WavePacket, ψₒ, ts)
	ts isa StepRangeLen || (ts = _toRange(ts))

	ψₒ = ψₒ isa Mode ? [ψₒ] : ψₒ

	length(sys) == 3 && return _ScatterProblem(sys[1], typeof(sys[1])[], sys[2], sys[3], ψᵢ, ψₒ, ts)
	return _ScatterProblem(sys[1], sys[2], sys[3], sys[4], ψᵢ, ψₒ, ts)
end
