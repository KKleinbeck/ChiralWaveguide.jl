mutable struct _ScatterProblem{O, S <: WavePacketState, WP <: _WavePacket{S}} <: WaveguideProblem
	H::O
	Ls::Array{O}
	σ::O
	system_state

	ψᵢ::WP
	ψₒ::Array{Mode}

	ts::StepRangeLen
end


function WaveguideProblem(sys, ψₒ, ts)
	WaveguideProblem(sys, NoPacket(), ψₒ, ts)
end

function WaveguideProblem(sys, ψᵢ::_WavePacket, ψₒ, ts)
	ts isa StepRangeLen || (ts = _toRange(ts))

	ψₒ = ψₒ isa Mode ? [ψₒ] : ψₒ

	length(sys) == 3 && return _ScatterProblem(sys[1], typeof(sys[1])[], sys[2], sys[3], ψᵢ, ψₒ, ts)
	return _ScatterProblem(sys[1], sys[2], sys[3], sys[4], ψᵢ, ψₒ, ts)
end
