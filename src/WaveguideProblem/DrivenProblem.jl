mutable struct _DrivenProblem{O, S <: WavePacketState, WP <: _WavePacket{S}} <: WaveguideProblem
	H::O
	Ls::Array{O}
	σ::O
	system_state

	ψᵢ::WP

	ts::StepRangeLen
end


function WaveguideProblem(sys, ts, disp_out::Bool = false)
	WaveguideProblem(sys, NoPacket(), ts, false)
end

function WaveguideProblem(sys, ψᵢ::_WavePacket, ts, disp_out::Bool = false)
	disp_out && @warn("Displacing the output ineffective without output modes")
	ts isa StepRangeLen || (ts = _toRange(ts))

	length(sys) == 3 && return _DrivenProblem(sys[1], typeof(sys[1])[], sys[2], sys[3], ψᵢ, ts)
	return _DrivenProblem(sys[1], sys[2], sys[3], sys[4], ψᵢ, ts)
end
