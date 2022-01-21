mutable struct _ScatterProblem{O, S <: WavePacketState, WP <: _WavePacket{S}} <: WaveguideProblem
	H::O
	Ls::Array{O}
	σ::O
	system_state

	ψᵢ::WP
	ψₒ::Array{Mode}

	ts::StepRangeLen

	displace_output::Bool
end


function WaveguideProblem(sys, ψₒ, ts, disp_out = false)
	WaveguideProblem(sys, NoPacket(), ψₒ, ts, disp_out)
end

function WaveguideProblem(sys, ψᵢ::_WavePacket, ψₒ, ts, disp_out = false)
	ts isa StepRangeLen || (ts = _toRange(ts))

	ψₒ = ψₒ isa Mode ? [ψₒ] : ψₒ

	if length(sys) == 3
		return _ScatterProblem(sys[1], typeof(sys[1])[], sys[2], sys[3], ψᵢ, ψₒ, ts, disp_out)
	end
	return _ScatterProblem(sys[1], sys[2], sys[3], sys[4], ψᵢ, ψₒ, ts, disp_out)
end
