"""
    DissipativeLambdaChain(N_atoms::Int; γd = 1.0, Γ = 0.0)

Creates a chain of `N_atoms` three-level atoms with spontaneous decay rate `Γ` and decay rate
`γd` into the second ground state.
The coupling rate is set to 1, but can be modified by scaling the Hamiltonian and coupling operator.
Uses `NLevelBasis(3)` as the single atom basis.
"""
function DissipativeLambdaChain(N_atoms::Int; γd = 0.0, Γ = 0.0)
	singleAtomBasis = NLevelBasis(3)
	σ_WG_i, σ_GW_i  = transition(singleAtomBasis, 2, 1), transition(singleAtomBasis, 1, 2)
	σ_DW_i = transition(singleAtomBasis, 3, 2)

	if N_atoms == 1
		groundstateAtoms = nlevelstate(singleAtomBasis, 1)
		H_Atom = 0.0 * one(singleAtomBasis)
		Dissipators = iszero(Γ) ? [√(γd) * σ_DW_i] : [√(γd) * σ_DW_i, √(Γ) * σ_GW_i]

		return H_Atom, Dissipators, σ_GW_i, groundstateAtoms
	end

	multiAtomBasis = reduce(⊗, repeat([singleAtomBasis], N_atoms))

	groundstateAtoms = reduce(⊗, repeat([nlevelstate(singleAtomBasis, 1)], N_atoms))

	σ_col = sum(i -> embed(multiAtomBasis, i, σ_GW_i), [1:N_atoms;])

	Dissipators = [embed(multiAtomBasis, i, √(γd) * σ_DW_i) for i ∈ 1:N_atoms]
	if !iszero(Γ)
		push!(Dissipators, [embed(multiAtomBasis, i, √(Γ) * σ_GW_i) for i ∈ 1:N_atoms]...)
	end

	H_Atoms = -0.5im * sum(i -> sum(
			j -> embed(multiAtomBasis, i, σ_WG_i) * embed(multiAtomBasis, j, σ_GW_i),
			[1:i-1;]
		), [2:N_atoms;]
	)
	H_Atoms = (H_Atoms + dagger(H_Atoms))

	return H_Atoms, Dissipators, σ_col, groundstateAtoms
end
