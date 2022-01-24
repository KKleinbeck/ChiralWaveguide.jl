# Displaced States

*ChiralWaveguide.jl* is designed with efficiency in mind.
One of the tricks used to make calculations more efficient is the Mollow transform.
The Mollow transform allows to take the displacement out of every displaced state (e.g., [`DisplacedFock`](@ref)) and put it into the Hamiltonian.
This drastically reduces the size of the input Hilbert space - for a coherent state we can even remove the input cavity entirely.
If for some reason, the full input cavity needs to be included, this can be manually achieved, at the loss of performance of course.

In this example we discuss [`ContinuousWave`](@ref), a special implementation of coherent states which generally provide the fastest implementation of a constant coherent state input.
We will compare this to alternative implementations, to discuss their particular implementation and use cases.

## Continuous wave input
Let's simulate a two level atom driven by a constant coherent light field.
We use [`ContinuousWave`](@ref) for this.
Using [`Coherent`](@ref) states works as well, and we will discuss this below.
It should be noted, however, that coherent states should always describe proper modes, i.e., normalisable wave functions, while continuous waves relax this condition.

To compare the performance with the other implementations we use `BenchmarkTools` `@btime` macro
```
tf = 5.0
α  = 1.0

problem = WaveguideProblem(TwoLevelChain(1), ContinuousWave(α), tf)
ts, ρs = @btime solve($problem) # 78.560 μs (367 allocations: 39.05 KiB)
```
`ContinuousWave` does two things to make the simulation fast:
Firstly, it tells the solver that the input is a coherent state and therefore the input cavity can be eliminated from the start.
This can be seen by
```
basis(ρs[end]) # NLevel(N=2)
```
i.e., `ρs` only stores the state of the atom.
Secondly, if a constant is passed to `ContinuousWave` the solver knows that the master equation is time independent and falls back to a faster implementation for time-independent master equations.
Of course `ContinuousWave` can also describe time dependence, but loosing the second advantage by doing so:
```
problem = WaveguideProblem(TwoLevelChain(1), ContinuousWave(t -> α), tf)
ts, ρs = @btime solve($problem) # 1.157 ms (20997 allocations: 1.84 MiB)
```
Already in this minimalistic example the difference is substantial.

## Coherent state input
Using `Coherent` instead of `ContinuousWave` primarily expresses the intent of using a proper mode, instead of a non-normalisable driving function.
Yet, the solver internally only uses the mode function of the specified mode, so it is possible to implement a continuous input with `Coherent`
```
problem = WaveguideProblem(
  TwoLevelChain(1),
  WavePacket(Mode(t -> 1., t -> 0., t -> 0., t -> 0.), Coherent(α)),
  tf
)
ts, ρs = @btime solve($problem) # 1.105 ms (20997 allocations: 1.84 MiB)
```
This provides the same performance as the time-independent implementation of `ContinuousWave`; in fact, this is how `ContinuousWave` is implemented for time dependent arguments.
However, this manual hack is discouraged, as `ContinuousWave` provides a simpler interface, that clearly expresses its intend and is potentially faster for constant driving amplitudes.
For normalisable wave functions both `Coherent` and `ContinuousWave` are practically interchangeable, but it is still recommended to use `Coherent`, as in the last example, to express clear intend.

## Accessing the input cavity
All previous examples eliminate the input cavity in favour for better performance.
If one wants to access the precise population of the input mode however, we can save the coherent number distribution in [`ArbitraryState`](@ref).
This needs a little boilerplating
```
N_cutoff = ceil(Int, tf * α^2 + 5 * √(tf) * α)
amplitudes = coherentstate(FockBasis(N_cutoff), √(tf) * α).data
state = ArbitraryState(amplitudes)
```
First, a cutoff is choosen, so that the coherent state fits reasonably well in the Hilbert space (here `sum(abs2, amplitudes) == 0.9999...`).
Then we use *QuantumOptics.jl* `coherentstate` to get the amplitudes and create our input state.
Notice, we choose an amplitude `√(tf) * α` and not `α`, since we use `HardBoxMode` as our mode function, which introduces a normalisation factor `1/√(tf)`.
Here it is important that we use a proper normalisable mode, since we need to access `gᵢ`.
The simulation now becomes
```
problem = WaveguideProblem(TwoLevelChain(1), WavePacket(HardBoxMode(σ = tf), state), tf)
ts, ρs_full = @btime solve($problem) # 207.809 ms (138290 allocations: 516.25 MiB)
```
and is significantly the slowest version.
Yet, we now may observe the input cavity and its entanglement with the quantum system or output cavity
```
basis(ρs_full[end]) # [Fock(cutoff=17) ⊗ NLevel(N=2)]
```
Lastly, we verify that this approach yields the correct result for the atom's state
```
ρs_atoms = [ptrace(ρ, 1) for ρ ∈ ρs_full]
all(tracedistance.(ρs_atoms, ρs) .< 2e-5) # true
```

## Mollow transformation on the output cavity
Let us now introduce an output cavity to the problem.
If we assume for the moment, that there is no quantum system between the input and output cavity,
then the time evolution of the state ``\rho`` of the output cavity follows the master equation
```math
\begin{aligned}
	\partial_t \rho   &= -i[H_\mathrm{drive}, \rho] + \mathcal{D}[\rho], \\
	H_\mathrm{drive}  &= i(\alpha^*(t) g_v^*(t) b - \alpha(t) g_v(t) b^\dagger), \\
	\mathcal{D}[\rho] &= |g_v(t)|^2 \Big[b \rho b^\dagger - \frac{1}{2} \{b^\dagger b, \rho\}\Big],
\end{aligned}
```
with ``g_v(t)`` the coupling rate of the output cavity.
The solution to this master equation for an initially de-excited cavity is a coherent state``\rho(t) = |\beta(t) \rangle\langle \beta(t)|`` with the amplitude
```math
\begin{equation}
  \beta(t) = - \int_0^t ds \alpha(s) g_v(s).
\end{equation}
```

This implies that, independent of the quantum system, we may expect that we need a large Fock basis for the output cavity of the order of ``|\alpha|^2``.
However, we may perform a Mollow transformation on the output cavity and eliminate ``H_\mathrm{drive}`` from the Master equation and thus save precious memory and make the simulations faster.
This becomes especially attractive when we are only interested in the final state ``\rho(t_f)``, when `\alpha(t_f) = 0` or `g_v(t_f) = 0`, so that we can exactly determine the displacement.
For example, for constant ``\alpha`` and a flat mode of width ``\sigma`` we find ``\beta(t_f) = \sqrt{\sigma} \alpha``.

Okay, let's see this trick in action.
In the simulation we can tell *ChiralWaveguide.jl* to do the Mollow transform by setting the `displace_output` flag to `true`:
```
probNonDisp = WaveguideProblem(TwoLevelChain(1), ContinuousWave(α), HardBoxMode(σ = tf), tf)
probDisp    = WaveguideProblem(TwoLevelChain(1), ContinuousWave(α), HardBoxMode(σ = tf), tf, true)

ts, ρsNonDisp = @btime solve($probNonDisp) # 121.945 ms (144023 allocations: 357.63 MiB)
ts, ρsDisp    = @btime solve($probDisp)    # 16.258 ms (57298 allocations: 53.28 MiB)
```
The strong improvement is of course due to the smaller Fock basis in the displaced case
```
ρCavNonDisp = ptrace(ρsNonDisp[end], 1)
ρCavDisp    = ptrace(ρsDisp[end],    1)

basis(ρCavNonDisp) # Fock(cutoff=14)
basis(ρCavDisp)    # Fock(cutoff=7)
```
Yet, If we account for the displacement the results coincide nicely
```
using GSL, SpecialFunctions
function displaceProper(basis, α)
	operator = DenseOperator(basis)

	αc, α2, eα2 = -conj(α), abs2(α), exp(-0.5abs2(α))

	for m ∈ 1:basis.N
		operator.data[m, m] = eα2 * sf_laguerre_n(m - 1, 0, α2)
		for n ∈ m:basis.N+1
			c = √(gamma(m)/gamma(n)) * eα2 * sf_laguerre_n(m - 1, n-m, α2)
			operator.data[m, n] = c * αc^(n-m)
			operator.data[n, m] = c * α^(n-m)
		end
	end
	return operator
end

D = displaceProper(basis(ρCavNonDisp), sqrt(tf) * α)
sqrt(sum(abs2, (D' * ρCavNonDisp * D).data[1:8,1:8] - ρCavDisp.data)) # 0.0023699962743152143
```
Manually changing `Nouts` in `solve` reveals that the difference is mainly due to the basis size of the *non displaced* implementation, i.e., the Mollow transformation is even more accurate.
Finally notice, we implemented our own version of the displacement operator and did not use `displace` from [QuantumOptics.jl](https://qojulia.org/), as this version do not provide correct matrix coefficient for any truncated Fock basis (it is however unitary, contrary to our implementation).
