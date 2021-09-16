"""
WaveguideProblem

Contains all parameters for the simulation of the chiral waveguide.
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
