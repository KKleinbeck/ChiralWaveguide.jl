mutable struct _DrivenProblem{O, S <: WavePacketState} <: WaveguideProblem
	H::O
	Ls::Array{O}
	σ::O
	system_state

	ψᵢ::WavePacket{S}

	ts::StepRangeLen
end

function WaveguideProblem(sys, ψᵢ::WavePacket, ts)
	ts isa StepRangeLen || (ts = _toRange(ts))

	length(sys) == 3 && return _DrivenProblem(sys[1], typeof(sys[1])[], sys[2], sys[3], ψᵢ, ts)
	return _DrivenProblem(sys[1], sys[2], sys[3], sys[4], ψᵢ, ts)
end
