# Defining the quantum system

Each quantum system used in a simulation consists of three ingredients:
1. The generator of the system's dynamic in the absence of the input cavity.
2. The operator that couples the external bosonic modes to the quantum system.
3. The initial state of the system.
These parameters are simply passed to `WaveguideProblem` as an array
```
WaveguideProblem([hamiltonian, σ, ψ₀], WavePacket(...), t)
```
or, for dissipative dynamic,
```
WaveguideProblem([hamiltonian, [Ls], σ, ψ₀], WavePacket(...), t)
```
It should be noted that, right now, the initial state of the quantum system
cannot be entangled with the input or output modes.

As a practical example we consider a two level atom, which we model as a spin-1/2 particle.
The atom has a dephasing `δ`, spontaneous decay rate `Γ`, starts in the excited state, and
is coupled to the bosonic modes via the interaction Hamiltonian
```math
H = - \frac{i}{2} \big[σ^+ b(0) - σ^- b^\dagger(0)].
```
The operator that couples the quantum system to the bosons is defined to be the operator next to
the ``b^\dagger(0)`` term, i.e., ``σ^-`` here.
To simulate this system we use [QuantumOptics.jl](https://qojulia.org/)'s spin operators
```
basis  = SpinBasis(1//2)
σ⁺, σ⁻ = sigmap(basis), sigmam(basis)

δ, Γ = rand(2)

H  = δ * σ⁺ * σ⁻
L  = √(Γ) * σ⁻
ψ₀ = spinup(basis)
```
Then, we simply create the `WaveguideProblem`
```
WaveguideProblem([H, [L], σ⁻, ψ₀], WavePacket(...), t)
```
By default *ChiralWaveguide.jl* assumes that the interaction Hamiltonian between bosons and the
quantum system has coupling rate ``1``.
To change this you should rescale all other rates (including the bosons' rates) by the respective
coupling rate you wish for.
Lastly, nothing is special about the spin-1/2 representation we choose here.
In fact the exported systems [`TwoLevelChain`](@ref) and [`DissipativeLambdaChain`](@ref) use
`NLevelBasis` instead and every basis in the *QuantumOptics.jl* ecosystem works.
