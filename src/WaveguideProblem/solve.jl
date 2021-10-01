"""
solve(problem::WaveguideProblem, args...; [Nouts], kwargs...)

Solve the `WaveguideProblem`. Additional arguments are passed to the time-dependant master equation
solver from the `QuantumOptics` package.
If there are output modes, then `NOuts` may specify the dimension of the respective Fockspaces.

Returns:
- ts:    The times on which the output is determined.
- ρ(t):  The density matrix at each point in time. The Hilbert space depends on the
         parameters given, if unsure about the result check out `ρ.basis_l`
"""
function solve(problem::WaveguideProblem; Nouts::Union{Nothing,Array} = nothing, kwargs...)
	if problem isa _ScatterProblem
		isnothing(Nouts) && ( Nouts = _expectedOutputPhotons(problem) )
		@assert length(Nouts) == length(problem.ψₒ)
	end

	ψ₀, Ls, H_nh = _generateLindbladian(problem, Nouts)

	function Operators(t, ψ)
		Ls_t = [Ls[1](t), Ls[2:end]...]
		return H_nh(t), H_nh(t)', Ls_t, dagger.(Ls_t)
	end

	return timeevolution.master_nh_dynamic(problem.ts, ψ₀, Operators; kwargs...)
end


function _expectedOutputPhotons(problem::_ScatterProblem{O, Coherent, ContinuousWave}) where {O}
	fill(coherent_cutoff( (problem.ts[end] - problem.ts[1]) * abs2(problem.ψᵢ.α) ),
		length(problem.ψₒ)
	)
	# For each mode, perform the integration to estimate the expectated number of photons
end
function _expectedOutputPhotons(problem::_ScatterProblem{O, S, WavePacket{S}}) where {O, S <: Displaced}
	N = problem.ψᵢ.state.N_cutoff + coherent_cutoff(problem.ψᵢ.state.α)
	fill(N, length(problem.ψₒ))
end
function _expectedOutputPhotons(problem::_ScatterProblem{O, S, WavePacket{S}}) where {O, S}
	N = problem.ψᵢ.state.N_cutoff
	fill(N, length(problem.ψₒ))
end


# ------------------------------------------------------------
# _DrivenProblems

function _generateLindbladian(problem::_DrivenProblem{O, Coherent, ContinuousWave}, Nouts) where {O}
	_generateLindbladian(
		WaveguideProblem(
			[problem.H, problem.Ls, problem.σ, problem.system_state],
			WavePacket(FlatMode(), Coherent(problem.ψᵢ.α)),
			problem.ts
		),
		Nouts
	)
end


function _generateLindbladian(problem::_DrivenProblem{O, Coherent, WavePacket{Coherent}}, Nouts) where {O}
	# ----------------------------------------
	# Scatter System
	H_sys, σ⁻, σ⁺ = problem.H, problem.σ, problem.σ'
	σ⁺σ⁻ = σ⁺ * σ⁻

	# ----------------------------------------
	# Constructing Disspator and Hamiltonian
	α  = problem.ψᵢ.state.α
	mf = problem.ψᵢ.mode.modeFunction
	atomic_dissipation = isempty(problem.Ls) ?
		-0.5im * σ⁺σ⁻ :
		-0.5im * σ⁺σ⁻ - 0.5im * sum([L'*L for L ∈ problem.Ls])

	H_nh(t) = atomic_dissipation - 1.0im * (σ⁺ * α * mf(t) - σ⁻ * α * mf(t)) + H_sys

	ψ₀ = problem.system_state
	return ψ₀, [t -> σ⁻, problem.Ls...], H_nh
end


function _generateLindbladian(problem::_DrivenProblem{O, T, WavePacket{T}}, Nouts) where {O, T}
	# ----------------------------------------
	# Scatter System
	H_sys, σ⁻, σ⁺ = problem.H, problem.σ, problem.σ'
	σ⁺σ⁻ = σ⁺ * σ⁻
	idₛ = one(H_sys.basis_l)

	# ----------------------------------------
	# Input System
	inBasis = FockBasis(problem.ψᵢ.state.N_cutoff)
	gᵢ      = problem.ψᵢ.mode.gᵢ
	idᵢ     = identityoperator(inBasis)
	aᵢ, aᵢᵀ = destroy(inBasis) ⊗ idₛ, create(inBasis) ⊗ idₛ
	aᵢᵀaᵢ   = aᵢᵀ * aᵢ

	# ----------------------------------------
	# Embed Operators and abbreviations
	H_sys, σ⁻, σ⁺, σ⁺σ⁻ = idᵢ ⊗ H_sys, idᵢ ⊗ σ⁻, idᵢ ⊗ σ⁺, idᵢ ⊗ σ⁺σ⁻
	aᵢᵀσ⁻,  aᵢσ⁺ = aᵢᵀ * σ⁻, aᵢ * σ⁺

	# ----------------------------------------
	# Constructing Disspator and Hamiltonian
	L  = t -> σ⁻ + gᵢ(t)aᵢ
	Lᵀ = t -> σ⁺ + gᵢ(t)aᵢᵀ

	atomic_dissipation = isempty(problem.Ls) ?
		-0.5im * σ⁺σ⁻ :
		-0.5im * σ⁺σ⁻ - (0.5im * idᵢ) ⊗ sum([L'*L for L ∈ problem.Ls])

	H_nh = _createHamiltonian(problem, atomic_dissipation, aᵢᵀaᵢ, aᵢσ⁺, σ⁺, σ⁻, H_sys)

	ψ₀ = Ket(inBasis, createState(problem.ψᵢ.state)) ⊗ problem.system_state
	return ψ₀, [L, [idᵢ ⊗ L for L ∈ problem.Ls]...], H_nh
end


function _createHamiltonian(problem::_DrivenProblem{O, S, WavePacket{S}},
		atomic_dissipation, aᵢᵀaᵢ, aᵢσ⁺, σ⁺, σ⁻, H_sys) where {O, S <: Displaced}
	gᵢ, mf, α = problem.ψᵢ.mode.gᵢ, problem.ψᵢ.mode.modeFunction, problem.ψᵢ.state.α

	return t -> atomic_dissipation - 0.5im * gᵢ(t)^2 * aᵢᵀaᵢ -
		1.0im * (σ⁺ * α * mf(t) - σ⁻ * α' * mf(t)') - 1.0im * gᵢ(t)aᵢσ⁺ + H_sys
end


function _createHamiltonian(problem::_DrivenProblem{O, S, WavePacket{S}},
		atomic_dissipation, aᵢᵀaᵢ, aᵢσ⁺, σ⁺, σ⁻, H_sys) where {O, S <: NonDisplaced}
	gᵢ = problem.ψᵢ.mode.gᵢ

	return t -> atomic_dissipation - 0.5im * gᵢ(t)^2 * aᵢᵀaᵢ - 1.0im * gᵢ(t)aᵢσ⁺ + H_sys
end


# ------------------------------------------------------------
# _ScatterProblems

function _generateLindbladian(problem::_ScatterProblem{O, Coherent, ContinuousWave}, Nouts) where {O}
	_generateLindbladian(
		WaveguideProblem(
			[problem.H, problem.Ls, problem.σ, problem.system_state],
			WavePacket(FlatMode(), Coherent(problem.ψᵢ.α)),
			problem.ψₒ,
			problem.ts
		),
		Nouts
	)
end

function _generateLindbladian(problem::_ScatterProblem{O, Coherent, WavePacket{Coherent}}, Nouts) where {O}
	# ----------------------------------------
	# Output System
	outBasis, groundstateOutput, aₒs, aₒᵀs = _generateOutputCavities(Nouts)
	idₒ = one(outBasis)

	# ----------------------------------------
	# Scatter System
	H_sys, σ⁻, σ⁺ = problem.H ⊗ idₒ, problem.σ ⊗ idₒ, problem.σ' ⊗ idₒ
	σ⁺σ⁻ = σ⁺ * σ⁻
	idₛ = one(problem.H.basis_l)

	# ----------------------------------------
	# Embed Operators and abbreviations
	aₒs, aₒᵀs = [idₛ ⊗ a for a ∈ aₒs], [idₛ ⊗ aᵀ for aᵀ ∈ aₒᵀs]

	gₒs = [mode.gₒ for mode in problem.ψₒ]

	gₒaₒs(t)  = sum([gₒ(t) for gₒ ∈ gₒs] .* aₒs)
	gₒaₒᵀs(t) = sum([gₒ(t) for gₒ ∈ gₒs] .* aₒᵀs)

	aₒᵀσ⁻s = [aₒᵀ * σ⁻ for aₒᵀ in aₒᵀs]
	aₒᵀaₒs = aₒᵀs .* aₒs
	# TODO: maybe it is faster to do `gₒaₒs(t)  = idₛ ⊗ sum(gₒs(t) .* aₒs)` and leave aₒs unchanged
	# due to the sum

	# ----------------------------------------
	# Constructing Disspator and Hamiltonian
	L  = t -> (σ⁻ + gₒaₒs(t) )
	Lᵀ = t -> (σ⁺ + gₒaₒᵀs(t) )

	α  = problem.ψᵢ.state.α
	mf = problem.ψᵢ.mode.modeFunction
	atomic_dissipation = isempty(problem.Ls) ?
		-0.5im * σ⁺σ⁻ :
		-0.5im * σ⁺σ⁻ - 0.5im * sum([L'*L for L ∈ problem.Ls]) ⊗ idₒ

	function H_nh(t)
		return atomic_dissipation - 0.5im * sum([gₒ(t)^2 for gₒ ∈ gₒs] .* aₒᵀaₒs) -
			1.0im * sum([gₒ(t) for gₒ ∈ gₒs] .* aₒᵀσ⁻s) -
			1.0im * (Lᵀ(t) * α * mf(t) - L(t) * α' * mf(t)) +
			H_sys
	end

	ψ₀ = problem.system_state ⊗ groundstateOutput
	return ψ₀, [t -> L(t), [L ⊗ idₒ for L ∈ problem.Ls]...], H_nh
end


function _generateLindbladian(problem::_ScatterProblem{O, T, WavePacket{T}}, Nouts) where {O, T}
	# ----------------------------------------
	# Output System
	outBasis, groundstateOutput, aₒs, aₒᵀs = _generateOutputCavities(Nouts)
	idₒ = one(outBasis)

	# ----------------------------------------
	# Scatter System
	H_sys, σ⁻, σ⁺  = problem.H, problem.σ, problem.σ'
	σ⁺σ⁻ = σ⁺ * σ⁻
	idₛ = one(problem.H.basis_l)

	# ----------------------------------------
	# Input System
	inBasis = FockBasis(problem.ψᵢ.state.N_cutoff)
	gᵢ      = problem.ψᵢ.mode.gᵢ
	idᵢ     = identityoperator(inBasis)
	aᵢ, aᵢᵀ = destroy(inBasis) ⊗ idₛ ⊗ idₒ, create(inBasis) ⊗ idₛ ⊗ idₒ
	aᵢᵀaᵢ = aᵢᵀ * aᵢ

	# ----------------------------------------
	# Embed Operators and abbreviations
	H_sys, σ⁻, σ⁺, σ⁺σ⁻ = idᵢ ⊗ H_sys ⊗ idₒ, idᵢ ⊗ σ⁻ ⊗ idₒ, idᵢ ⊗ σ⁺ ⊗ idₒ, idᵢ ⊗ σ⁺σ⁻ ⊗ idₒ
	aₒs, aₒᵀs = [idᵢ ⊗ idₛ ⊗ a for a ∈ aₒs], [idᵢ ⊗ idₛ ⊗ aᵀ for aᵀ ∈ aₒᵀs]

	gₒs = [mode.gₒ for mode in problem.ψₒ]

	gₒaₒs(t)  = sum([gₒ(t)' for gₒ ∈ gₒs] .* aₒs)
	gₒaₒᵀs(t) = sum([gₒ(t)  for gₒ ∈ gₒs] .* aₒᵀs)

	aᵢσ⁺   = aᵢ * σ⁺
	aₒᵀaᵢs = [aₒᵀ * aᵢ for aₒᵀ in aₒᵀs]
	aₒᵀσ⁻s = [aₒᵀ * σ⁻ for aₒᵀ in aₒᵀs]
	aₒᵀaₒs = aₒᵀs .* aₒs
	# TODO: maybe it is faster to do `gₒaₒs(t)  = idₛ ⊗ sum(gₒs(t) .* aₒs)` and leave aₒs unchanged
	# due to the sum

	# ----------------------------------------
	# Constructing Disspator and Hamiltonian
	L  = t -> σ⁻ + gᵢ(t)aᵢ  + gₒaₒs(t)

	atomic_dissipation = isempty(problem.Ls) ?
		-0.5im * σ⁺σ⁻ :
		-0.5im * σ⁺σ⁻ - (0.5im * idᵢ) ⊗ sum([L'*L for L ∈ problem.Ls]) ⊗ idₒ

	H_nh = _createHamiltonian(problem, atomic_dissipation,
		aₒᵀaₒs, aᵢᵀaᵢ, aᵢσ⁺, aₒᵀaᵢs, aₒᵀσ⁻s, σ⁺, σ⁻, gₒaₒᵀs, gₒaₒs, H_sys
	)

	ψ₀ = Ket(inBasis, createState(problem.ψᵢ.state)) ⊗ problem.system_state ⊗ groundstateOutput
	return ψ₀, [t -> L(t), [idᵢ ⊗ L ⊗ idₒ for L ∈ problem.Ls]...], H_nh
end


function _generateOutputCavities(Nouts)
	# Preconverting into `CompositeBasis` allows us to `embed` into it, even for a single system
	outBasis  = reduce(⊗, FockBasis.(Nouts) .|> CompositeBasis)

	groundstateOutput = reduce(⊗, [fockstate(FockBasis(N), 0) for N in Nouts])

	aₒs  = [embed(outBasis, i, destroy(FockBasis(N))) for (i, N) in enumerate(Nouts)]
	aₒᵀs = [embed(outBasis, i, create(FockBasis(N)))  for (i, N) in enumerate(Nouts)]

	return outBasis, groundstateOutput, aₒs, aₒᵀs
end


function _createHamiltonian(problem::_ScatterProblem{O, S, WavePacket{S}},
		atomic_dissipation, aₒᵀaₒs, aᵢᵀaᵢ, aᵢσ⁺, aₒᵀaᵢs, aₒᵀσ⁻s, σ⁺, σ⁻, gₒaₒᵀs, gₒaₒs,
		H_sys) where {O, S<: Displaced}
	gᵢ, mf, α = problem.ψᵢ.mode.gᵢ, problem.ψᵢ.mode.modeFunction, problem.ψᵢ.state.α
	gₒs   = [mode.gₒ for mode in problem.ψₒ]

	return t -> atomic_dissipation -
		0.5im * (gᵢ(t)^2 * aᵢᵀaᵢ + sum([gₒ(t)^2 for gₒ ∈ gₒs] .* aₒᵀaₒs) ) -
		1.0im * (gᵢ(t)aᵢσ⁺ + sum([gₒ(t) for gₒ ∈ gₒs] .* aₒᵀσ⁻s) +
			sum(gᵢ(t) .* [gₒ(t) for gₒ ∈ gₒs] .* aₒᵀaᵢs) ) -
		1.0im * (( σ⁺ + gₒaₒᵀs(t) ) * α * mf(t) - ( σ⁻ + gₒaₒs(t) )  * α' * mf(t)') + H_sys
end


function _createHamiltonian(problem::_ScatterProblem{O, S, WavePacket{S}},
		atomic_dissipation, aₒᵀaₒs, aᵢᵀaᵢ, aᵢσ⁺, aₒᵀaᵢs, aₒᵀσ⁻s, σ⁺, σ⁻, gₒaₒᵀs, gₒaₒs,
		H_sys) where {O, S<: NonDisplaced}
	gᵢ    = problem.ψᵢ.mode.gᵢ
	gₒs   = [mode.gₒ for mode in problem.ψₒ]

	return t -> atomic_dissipation -
		0.5im * (gᵢ(t)^2 * aᵢᵀaᵢ + sum([gₒ(t)^2 for gₒ ∈ gₒs] .* aₒᵀaₒs) ) -
		1.0im * (gᵢ(t)aᵢσ⁺ + sum([gₒ(t) for gₒ ∈ gₒs] .* aₒᵀσ⁻s) +
			sum(gᵢ(t) .* [gₒ(t) for gₒ ∈ gₒs] .* aₒᵀaᵢs) ) +
		H_sys
end
