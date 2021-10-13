"""
    WaveguideProblem(system, [ψᵢ::WavePacket, ψ₀::Mode,] ts)

Defines the waveguide problem for [`solve(x)`](@ref).

# Arguments
- system: array of form `[H, σ, ψₛ]` or `[H, Ls, σ, ψₛ]`. `H` is the system Hamiltonian, `Ls` an array
  of system dissipators, `σ` the coupling operator of the system, and `ψₛ` the systems initial state.
- ψᵢ (optional): the [`WavePacket`](@ref) driving the system.
- ψₒ (optional): the observed mode.
- ts: time domain for the simulation. Can be a number, tuple or array. If `ts` is a number the
  simulation is started at `t = 0`.
"""
abstract type WaveguideProblem end

function _toRange(t::Number)
	@assert t > 0 "Default initial time t₀ = 0 does not work with negative final times"
	range(0, t, length = 101)
end
function _toRange(ts::Tuple{N,N}) where N <: Number
	@assert ts[1] < ts[2] "Times are not ordered"
	range(ts[1], ts[2], length = 101)
end

include("./DrivenProblem.jl")

include("./ScatterProblem.jl")
