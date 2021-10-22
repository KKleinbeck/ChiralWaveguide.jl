# Simulating 1d chiral quantum systems

*ChiralWaveguide.jl* is a simulation package for arbitrary chiral quantum systems, coupled to a
continuous bosonic background with linear dispersion, for example photons.
This package aims to deliver three promises:

1. A lightweight, descriptive interface. Scattering a Gaussian Fock state on a chiral atom is no
   more complicated than
   ```
   problem = WaveguideProblem(TwoLevelChain(1), WavePacket(GaussMode(), Fock(1)), 10.0)
   ts, ρs = solve(problem)
   ```

2. Performance. *ChiralWaveguide.jl* performs suitable transformations on the Hilbert space, like
   the Mollow transform, to reduce the effective Hilbert space size as much as possible.

3. Interoperability and a native *Julia* experience. *ChiralWaveguide* is build entirely around
   the [QuantumOptics.jl](https://qojulia.org/) package. In fact, getting the population of the
   excited state in the previous example is no more complicated than
   ```
   σ⁺σ⁻ = transition(NLevelBasis(2), 2, 2)
   expect(2, σ⁺σ⁻, ρs)
   ```
   Additionally, the syntax is heavily inspired by the [SciML](https://sciml.ai/) ecosystem.
