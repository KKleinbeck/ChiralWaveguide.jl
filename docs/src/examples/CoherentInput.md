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
ts, ρs = @btime solve($problem) # 161.369 ms (21107 allocations: 2.88 MiB)
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
ts, ρs = @btime solve($problem) # 204.156 ms (51074 allocations: 4.77 MiB)
```
In the minimalistic example the difference is small, but becomes more significant for larger systems.

## Coherent state input
Using `Coherent` instead of `ContinuousWave` primarily expresses the intent of using a proper mode, instead of a non-normalisable driving function.
Yet, the solver internally only uses the mode function of the specified mode, so it is possible to implement a continuous input with `Coherent`
```
problem = WaveguideProblem(
  TwoLevelChain(1),
  WavePacket(Mode(t -> 1., t -> 0., t -> 0., t -> 0.), Coherent(α)),
  tf
)
ts, ρs = @btime solve($problem) # 210.791 ms (51074 allocations: 4.77 MiB)
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
ts, ρs_full = @btime solve($problem) # 1.300 s (282065 allocations: 560.71 MiB)
```
and is the significant slowest version.
Yet, we now may observe the input cavity and its entanglement with the quantum system or output cavity
```
basis(ρs_full[end]) # [Fock(cutoff=17) ⊗ NLevel(N=2)]
```
Lastly, we verify that this approach yields the correct result for the atom's state
```
ρs_atoms = [ptrace(ρ, 1) for ρ ∈ ρs_full]
all(tracedistance.(ρs_atoms, ρs) .< 2e-5) # true
```
