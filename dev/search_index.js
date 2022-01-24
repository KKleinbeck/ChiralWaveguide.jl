var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/#Solve","page":"API","title":"Solve","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"ChiralWaveguide.solve(x)","category":"page"},{"location":"api/#CommonSolve.solve-Tuple{Any}","page":"API","title":"CommonSolve.solve","text":"solve(problem::WaveguideProblem; [ρ₀, Nouts], kwargs...)\n\nSolves the WaveguideProblem, using QuantumOptics.timeevolution.master_nh[_dynamic].\n\nArguments\n\nproblem: the waveguide problem\nNouts::Array{Int} (optional): specifies the dimension of the output Fock spaces\nρ₀ (optional): Alternative initial state for the simulation.\n\nOverrides the initial state provided by the `WavePacket` and the quantum system.\n\nAdditional keyword arguments are passed to QuantumOptics.jl's master equation solver, see it's documentation\n\nReturns\n\nts:    The times on which the output is determined.\nρ(t):  The density matrix at each point in time. The Hilbert space depends on the specific        problem; if unsure check out basis(ρ[1]).\n\n\n\n\n\n","category":"method"},{"location":"api/#WaveguideProblem","page":"API","title":"WaveguideProblem","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"ChiralWaveguide.WaveguideProblem","category":"page"},{"location":"api/#ChiralWaveguide.WaveguideProblem","page":"API","title":"ChiralWaveguide.WaveguideProblem","text":"WaveguideProblem(system, [ψᵢ::WavePacket, ψ₀::Array{Mode},] ts, [displace_output])\n\nDefines the waveguide problem for solve(x).\n\nArguments\n\nsystem: array of form [H, σ, ψₛ] or [H, Ls, σ, ψₛ]. H is the system Hamiltonian, Ls an array of system dissipators, σ the coupling operator of the system, and ψₛ the systems initial state.\nψᵢ (optional): the WavePacket driving the system.\nψₒ (optional): the observed modes (only one possible at the moment).\nts: time domain for the simulation. Can be a number, tuple or array. If ts is a number the simulation is started at t = 0.\ndisplace_output (optional): Removes the couplong between input and output modes when the input is a coherent state, i.e., work in a frame where the output modes are displaced accordingly.\n\n\n\n\n\n","category":"type"},{"location":"api/#WavePackets-and-States","page":"API","title":"WavePackets & States","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"WavePacket","category":"page"},{"location":"api/#ChiralWaveguide.WavePacket","page":"API","title":"ChiralWaveguide.WavePacket","text":"WavePacket(mode::Mode, state::WavePacketState)\n\nDescribes the input wave packet. A WavePacket consists of any Mode and any WavePacketState.\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"API","title":"API","text":"WavePacketState","category":"page"},{"location":"api/#ChiralWaveguide.WavePacketState","page":"API","title":"ChiralWaveguide.WavePacketState","text":"WavePacketState\n\nAbstract base class for all states.\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"API","title":"API","text":"Coherent","category":"page"},{"location":"api/#ChiralWaveguide.Coherent","page":"API","title":"ChiralWaveguide.Coherent","text":"Coherent(α[, N_cutoff])\n\nDescribes a coherent state α. The autmatic choice for N_cutoff yields a normalisation of above 0.999 for α ∈ [0,10].\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"API","title":"API","text":"ContinuousWave","category":"page"},{"location":"api/#ChiralWaveguide.ContinuousWave","page":"API","title":"ChiralWaveguide.ContinuousWave","text":"ContinuousWave(α::Union{Number, Function})\n\nA special form of a Coherent state. For constant α this allows to solver to perform additional optimisations. Additionally, expresses the intend of using a non-normalisable input mode.\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"API","title":"API","text":"ArbitraryState","category":"page"},{"location":"api/#ChiralWaveguide.ArbitraryState","page":"API","title":"ChiralWaveguide.ArbitraryState","text":"ArbitraryState(amplitudes[, N_cutoff])\n\nCreates a state with specific amplitudes.\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"API","title":"API","text":"DisplacedArbitraryState","category":"page"},{"location":"api/#ChiralWaveguide.DisplacedArbitraryState","page":"API","title":"ChiralWaveguide.DisplacedArbitraryState","text":"DisplacedArbitraryState(α, amplitudes[, N_cutoff])\n\nDescribes a state with specific amplitudes, which then is displaced by D(α).\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"API","title":"API","text":"Fock","category":"page"},{"location":"api/#ChiralWaveguide.Fock","page":"API","title":"ChiralWaveguide.Fock","text":"Fock(n[, N_cutoff])\n\nDecribes the Fock state n.\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"API","title":"API","text":"DisplacedFock","category":"page"},{"location":"api/#ChiralWaveguide.DisplacedFock","page":"API","title":"ChiralWaveguide.DisplacedFock","text":"  DisplacedFock(α, n[, N_cutoff])\n\nDecribes the displaced Fock state D(α)n.\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"API","title":"API","text":"SqueezedVacuum","category":"page"},{"location":"api/#ChiralWaveguide.SqueezedVacuum","page":"API","title":"ChiralWaveguide.SqueezedVacuum","text":"SqueezedVacuum(r, ϕ = 0.0[, N_cutoff])\n\tSqueezedVacuum(ξ::Complex{Float64}[, N_cutoff])\n\nDecribes the a squeezed state with squeezing amplitude r and squeezing angle ϕ, i.e., the state exp(ξ^* a^2 - ξ a^dagger 2)2 0 with ξ = r e^i ϕ.\n\n\n\n\n\n","category":"type"},{"location":"api/#Modes","page":"API","title":"Modes","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Mode","category":"page"},{"location":"api/#ChiralWaveguide.Mode","page":"API","title":"ChiralWaveguide.Mode","text":"Mode(modeFunction; compression = :algebraic, kwargs...)\nMode(modeFunction, gᵢ, gₒ, norm)\n\nContainer for a cavity mode.\n\nArguments\n\nmodeFunction: the wave function of the mode\ngᵢ: coupling function for an input cavity\ngₒ: coupling function for an output cavity\nnorm: square integral of the mode up to time t\n\nIf only the mode function is provided, then all other quantities are numerically determined. This is done by mapping ℝ  (-1 1) and then numerically solving the resulting integrals. The argument compression picks the map ℝ  (-1 1) and can be a Compression or the respective symbol. Allowed symbols are the keys of the ChiralWaveguide.Compressions dictionary. Additional keyword arguments will be passed to the numerical integrator.\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"API","title":"API","text":"FourierMode","category":"page"},{"location":"api/#ChiralWaveguide.FourierMode","page":"API","title":"ChiralWaveguide.FourierMode","text":"FourierMode(a0, an = [], bn = []; t₀ = 0.0, σ = 1.0)\n\nGeneral Fourier Series for a mode in the time bin (t₀, t₀ + σ), with mode function\n\nm(t) = a0  σ + frac2σ sum_n a_n cos(2π n (t - t₀)  σ) + frac2σ sum_n b_n sin(2π n (t - t₀)  σ)\n\nThe coefficients a0, an, bn do not have to be normalised to provide the correct coupling rates.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"API","title":"API","text":"HardBoxMode","category":"page"},{"location":"api/#ChiralWaveguide.HardBoxMode","page":"API","title":"ChiralWaveguide.HardBoxMode","text":"HardBoxMode(; t₀ = 0.0, σ = 1.0)\n\nThe mode is a box function, i.e., it takes the value 1  sqrtσ for t₀  t  t₀ + σ and 0 otherwise. This can be regarded as the special case of a FourierMode, i.e., FourierMode(1, [], []).\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"API","title":"API","text":"SoftBoxMode","category":"page"},{"location":"api/#ChiralWaveguide.SoftBoxMode","page":"API","title":"ChiralWaveguide.SoftBoxMode","text":"SoftBoxMode(; τ = 0.0, σ = 1.0, n::Int = 10)\n\nThe mode is an approximation to the box mode, mathrmbox(t)  frac1sqrt(2t)^2n + 1. The parameters τ and σ control the center and width respectively, n the exponent.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"API","title":"API","text":"SoftBoxExpMode","category":"page"},{"location":"api/#ChiralWaveguide.SoftBoxExpMode","page":"API","title":"ChiralWaveguide.SoftBoxExpMode","text":"SoftBoxExpMode(; τ = 0.0, σ = 1.0, γ = 10.0)\n\nA box mode with exponential decaying flanks, given by\n\n                   ⌜ 1                     for |t| < 1/2\nsoftBoxExp(t, γ) = |\n                   ⌞ exp(-γ(|t| - 1/2))    for |t| > 1/2\n\nThe parameters τ and σ control the center and width respectively, γ the deay rate.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"API","title":"API","text":"GaussMode","category":"page"},{"location":"api/#ChiralWaveguide.GaussMode","page":"API","title":"ChiralWaveguide.GaussMode","text":"GaussMode(; τ = 0.0, σ = 1.0)\n\nNumerical stable mode function & couplings of Gaussian wave packet with mean τ and variance σ.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"API","title":"API","text":"ExpMode","category":"page"},{"location":"api/#ChiralWaveguide.ExpMode","page":"API","title":"ChiralWaveguide.ExpMode","text":"ExpMode(; t₀ = 0.0, γ = 1.0)\n\nExponentially raises up to time t₀ or decays after t₀, depending on the sign of the rate γ. Amplitude decays with rate γ/2, so that the probability density decays with γ.\n\n\n\n\n\n","category":"function"},{"location":"api/#Systems","page":"API","title":"Systems","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"TwoLevelChain","category":"page"},{"location":"api/#ChiralWaveguide.TwoLevelChain","page":"API","title":"ChiralWaveguide.TwoLevelChain","text":"TwoLevelChain(N_atoms::Int; Γ = 0.0)\n\nCreates a chain of N_atoms two-level atoms with spontaneous decay rate Γ. The coupling rate is set to 1, but can be modified by scaling the Hamiltonian and coupling operator. Uses NLevelBasis(2) as the single atom basis.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"API","title":"API","text":"DissipativeLambdaChain","category":"page"},{"location":"api/#ChiralWaveguide.DissipativeLambdaChain","page":"API","title":"ChiralWaveguide.DissipativeLambdaChain","text":"DissipativeLambdaChain(N_atoms::Int; γd = 1.0, Γ = 0.0)\n\nCreates a chain of N_atoms three-level atoms with spontaneous decay rate Γ and decay rate γd into the second ground state. The coupling rate is set to 1, but can be modified by scaling the Hamiltonian and coupling operator. Uses NLevelBasis(3) as the single atom basis.\n\n\n\n\n\n","category":"function"},{"location":"api/#Compressions","page":"API","title":"Compressions","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"ChiralWaveguide.Compression","category":"page"},{"location":"api/#ChiralWaveguide.Compression","page":"API","title":"ChiralWaveguide.Compression","text":"Compression(c, d, i; σ = 1.0, μ = 0.0)\n\nMaps ℝ  (-1 1).\n\nArguments\n\nc: the \"compression\", i.e., the map ℝ  (-1 1),\nd: the \"decompression\", i.e., the map (-1 1)  ℝ,\ni: the jacobian of the map, i.e., d(τ) τ,\nσ, μ: control the width and offset in the numerical integrations\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"API","title":"API","text":"Algebraic","category":"page"},{"location":"api/#ChiralWaveguide.Algebraic","page":"API","title":"ChiralWaveguide.Algebraic","text":"Algebraic(; σ = 1.0, μ = 0.0)\n\nMaps ℝ  (-1 1) by t  fract1 + t. Equivalent symbol is :algebraic.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"API","title":"API","text":"Exponential","category":"page"},{"location":"api/#ChiralWaveguide.Exponential","page":"API","title":"ChiralWaveguide.Exponential","text":"Exponential(; σ = 1.0, μ = 0.0)\n\nMaps ℝ  (-1 1) by t  mathrmtanh(t). Equivalent symbol is :exponential.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"API","title":"API","text":"Trigonometric","category":"page"},{"location":"api/#ChiralWaveguide.Trigonometric","page":"API","title":"ChiralWaveguide.Trigonometric","text":"Trigonometric(; σ = 1.0, μ = 0.0)\n\nMaps ℝ  (-1 1) by t  frac2π mathrmtan^-1(t). Equivalent symbol is :trigonometric.\n\n\n\n\n\n","category":"function"},{"location":"theory/#Theory","page":"Theory","title":"Theory","text":"","category":"section"},{"location":"theory/","page":"Theory","title":"Theory","text":"To be filled.","category":"page"},{"location":"examples/CustomModes/#Creating-custom-modes","page":"Creating custom modes","title":"Creating custom modes","text":"","category":"section"},{"location":"examples/CustomModes/","page":"Creating custom modes","title":"Creating custom modes","text":"The constructor of Mode automatically determines the coupling rates","category":"page"},{"location":"examples/CustomModes/","page":"Creating custom modes","title":"Creating custom modes","text":"beginaligned\n\tgᵢ =  fracu(t)sqrt1 - int_-infty^t u(s)^2 ds \n\tgₒ = -fracu(t)sqrt    int_-infty^t u(s)^2 ds\nendaligned","category":"page"},{"location":"examples/CustomModes/","page":"Creating custom modes","title":"Creating custom modes","text":"from the provided wave function u, by numerically integration the integral in the denominator, stored in Mode.norm. This process is designed to be as robust as possible, however, for special choices of u the integrator may struggle. This is for example the case when u sharply peaked or peaked far away from the origin. In these situations we have two options:","category":"page"},{"location":"examples/CustomModes/","page":"Creating custom modes","title":"Creating custom modes","text":"Give the integrator hints about the wave function.\nProvide the analytic coupling rates.","category":"page"},{"location":"examples/CustomModes/","page":"Creating custom modes","title":"Creating custom modes","text":"We will now explore both options for the wave function u(t) = sqrtγexp(-γt-τ).","category":"page"},{"location":"examples/CustomModes/#Aiding-the-integrator","page":"Creating custom modes","title":"Aiding the integrator","text":"","category":"section"},{"location":"examples/CustomModes/","page":"Creating custom modes","title":"Creating custom modes","text":"ChiralWaveguide provides the tools to easily combat the most typical problems. If we set τ  1 in the upper example the integrator fails to find suitable integration points","category":"page"},{"location":"examples/CustomModes/","page":"Creating custom modes","title":"Creating custom modes","text":"u(t, τ, γ) = √(γ) * exp(-γ * abs(t - τ))\n\nτ = 0.0\npA = plot(τ .+ [-5:0.1:5;], Mode(t -> u(t, τ, 0)).norm, title = \"τ = $(τ)\", label = nothing)\n\nτ = 10.0\npB = plot(τ .+ [-5:0.1:5;], Mode(t -> u(t, τ, 0)).norm, title = \"τ = $(τ)\", label = nothing)\n\nplot(pA, pB, size = (900, 300))","category":"page"},{"location":"examples/CustomModes/","page":"Creating custom modes","title":"Creating custom modes","text":"(Image: Populations) To find Mode.norm the integrator first maps ℝ  (-1 1). The map can be specified by passing a ChiralWaveguide.Compression to the compression keyword argument. Here, we can also specify the approximate center and width of the wave function. The default choice is Algebraic, so lets try this here","category":"page"},{"location":"examples/CustomModes/","page":"Creating custom modes","title":"Creating custom modes","text":"plot(τ .+ [-5:0.1:5;], Mode(t -> u(t, τ, 0), compression = Algebraic(μ = τ)).norm,\n\ttitle = \"τ = $(τ), Aided\", label = nothing\n)","category":"page"},{"location":"examples/CustomModes/","page":"Creating custom modes","title":"Creating custom modes","text":"(Image: Populations)","category":"page"},{"location":"examples/CustomModes/","page":"Creating custom modes","title":"Creating custom modes","text":"Additionally, a sharply peaked wave function may decrease the accuracy of the determined norm and, in turn, of gᵢ and gₒ. Consider","category":"page"},{"location":"examples/CustomModes/","page":"Creating custom modes","title":"Creating custom modes","text":"τ = 0.0\nγ = 50\nMode(t -> u(t, τ, γ)).norm(1e5)\n# == 1.016513761511986","category":"page"},{"location":"examples/CustomModes/","page":"Creating custom modes","title":"Creating custom modes","text":"Even though u is normalised, the integrator over estimates the result. While the construction of gᵢ and gₒ are build with non-normalised wave functions in mind and are resilient against deviations (they won't diverge) this may still lead to small inaccuracies. In problems where accuracy is key, we want to tell the integrator the approximate width of the wave function","category":"page"},{"location":"examples/CustomModes/","page":"Creating custom modes","title":"Creating custom modes","text":"Mode(t -> u(t, τ, γ), compression = Algebraic(μ = τ, σ = 1/γ)).norm(1e5)\n# == 0.9990907077615423","category":"page"},{"location":"examples/CustomModes/","page":"Creating custom modes","title":"Creating custom modes","text":"This already yields 3 significant digits and likely is good enough for all applications. However, we can still improve upon this by increasing the accuracy of the integrator","category":"page"},{"location":"examples/CustomModes/","page":"Creating custom modes","title":"Creating custom modes","text":"Mode(t -> wf(t, τ, γ), compression = Algebraic(μ = τ, σ = 1/γ), reltol = 1e-6).norm(1e5)\n# == 0.999999763224909","category":"page"},{"location":"examples/CustomModes/","page":"Creating custom modes","title":"Creating custom modes","text":"Lastly, it should be noted that potential problems are manyfold and before using any of the described techniques one should always screen for the precise origin of the problem first. Especially mode functions which rely on special functions like the error function can be problematic since their implementation might not be stable at large t and one should probably introduce a manual cutoff.","category":"page"},{"location":"examples/CustomModes/#Providing-exact-coupling-rates","page":"Creating custom modes","title":"Providing exact coupling rates","text":"","category":"section"},{"location":"examples/CustomModes/","page":"Creating custom modes","title":"Creating custom modes","text":"Provided we have an analytic formula for int_-infty^t u(s)^2 ds we may also provide the rates gᵢ and gₒ directly. in our example we have","category":"page"},{"location":"examples/CustomModes/","page":"Creating custom modes","title":"Creating custom modes","text":"int_-infty^t u(s)^2 ds =\nbegincases\n\t    frac12exp(-2γt-τ)  t  τ \n\t1 - frac12exp(-2γt-τ)  t  τ\nendcases","category":"page"},{"location":"examples/CustomModes/","page":"Creating custom modes","title":"Creating custom modes","text":"The creation of the exact mode function is now accomplished with","category":"page"},{"location":"examples/CustomModes/","page":"Creating custom modes","title":"Creating custom modes","text":"function norm(t, τ, γ)\n\tt < τ && return 0.5 * exp(-2γ * (τ - t))\n\treturn 1 - 0.5 * exp(-2γ * (t - τ))\nend\n\ngᵢ(t, τ, γ) =  u(t, τ, γ) / √(1 - norm(t, τ, γ))\ngₒ(t, τ, γ) = -u(t, τ, γ) / √(    norm(t, τ, γ))\n\nτ, γ = 0.0, 1.0\nmode = Mode(t -> u(t, τ, γ), t -> gᵢ(t, τ, γ), t -> gₒ(t, τ, γ), t -> norm(t, τ, γ))","category":"page"},{"location":"examples/CustomModes/","page":"Creating custom modes","title":"Creating custom modes","text":"Note that if you only ever want to use your mode as a input (output) mode, it suffices to simply provide u and/or gᵢ (gₒ). norm is only used when the mode function is constructor automatically and can be ignored most of the time. The drawbacks of course are that you loose flexibility and possibly a way to do consistency checks, so use this power with care.","category":"page"},{"location":"examples/CoherentInput/#Displaced-States","page":"Displaced States","title":"Displaced States","text":"","category":"section"},{"location":"examples/CoherentInput/","page":"Displaced States","title":"Displaced States","text":"ChiralWaveguide.jl is designed with efficiency in mind. One of the tricks used to make calculations more efficient is the Mollow transform. The Mollow transform allows to take the displacement out of every displaced state (e.g., DisplacedFock) and put it into the Hamiltonian. This drastically reduces the size of the input Hilbert space - for a coherent state we can even remove the input cavity entirely. If for some reason, the full input cavity needs to be included, this can be manually achieved, at the loss of performance of course.","category":"page"},{"location":"examples/CoherentInput/","page":"Displaced States","title":"Displaced States","text":"In this example we discuss ContinuousWave, a special implementation of coherent states which generally provide the fastest implementation of a constant coherent state input. We will compare this to alternative implementations, to discuss their particular implementation and use cases.","category":"page"},{"location":"examples/CoherentInput/#Continuous-wave-input","page":"Displaced States","title":"Continuous wave input","text":"","category":"section"},{"location":"examples/CoherentInput/","page":"Displaced States","title":"Displaced States","text":"Let's simulate a two level atom driven by a constant coherent light field. We use ContinuousWave for this. Using Coherent states works as well, and we will discuss this below. It should be noted, however, that coherent states should always describe proper modes, i.e., normalisable wave functions, while continuous waves relax this condition.","category":"page"},{"location":"examples/CoherentInput/","page":"Displaced States","title":"Displaced States","text":"To compare the performance with the other implementations we use BenchmarkTools @btime macro","category":"page"},{"location":"examples/CoherentInput/","page":"Displaced States","title":"Displaced States","text":"tf = 5.0\nα  = 1.0\n\nproblem = WaveguideProblem(TwoLevelChain(1), ContinuousWave(α), tf)\nts, ρs = @btime solve($problem) # 161.369 ms (21107 allocations: 2.88 MiB)","category":"page"},{"location":"examples/CoherentInput/","page":"Displaced States","title":"Displaced States","text":"ContinuousWave does two things to make the simulation fast: Firstly, it tells the solver that the input is a coherent state and therefore the input cavity can be eliminated from the start. This can be seen by","category":"page"},{"location":"examples/CoherentInput/","page":"Displaced States","title":"Displaced States","text":"basis(ρs[end]) # NLevel(N=2)","category":"page"},{"location":"examples/CoherentInput/","page":"Displaced States","title":"Displaced States","text":"i.e., ρs only stores the state of the atom. Secondly, if a constant is passed to ContinuousWave the solver knows that the master equation is time independent and falls back to a faster implementation for time-independent master equations. Of course ContinuousWave can also describe time dependence, but loosing the second advantage by doing so:","category":"page"},{"location":"examples/CoherentInput/","page":"Displaced States","title":"Displaced States","text":"problem = WaveguideProblem(TwoLevelChain(1), ContinuousWave(t -> α), tf)\nts, ρs = @btime solve($problem) # 204.156 ms (51074 allocations: 4.77 MiB)","category":"page"},{"location":"examples/CoherentInput/","page":"Displaced States","title":"Displaced States","text":"In the minimalistic example the difference is small, but becomes more significant for larger systems.","category":"page"},{"location":"examples/CoherentInput/#Coherent-state-input","page":"Displaced States","title":"Coherent state input","text":"","category":"section"},{"location":"examples/CoherentInput/","page":"Displaced States","title":"Displaced States","text":"Using Coherent instead of ContinuousWave primarily expresses the intent of using a proper mode, instead of a non-normalisable driving function. Yet, the solver internally only uses the mode function of the specified mode, so it is possible to implement a continuous input with Coherent","category":"page"},{"location":"examples/CoherentInput/","page":"Displaced States","title":"Displaced States","text":"problem = WaveguideProblem(\n  TwoLevelChain(1),\n  WavePacket(Mode(t -> 1., t -> 0., t -> 0., t -> 0.), Coherent(α)),\n  tf\n)\nts, ρs = @btime solve($problem) # 210.791 ms (51074 allocations: 4.77 MiB)","category":"page"},{"location":"examples/CoherentInput/","page":"Displaced States","title":"Displaced States","text":"This provides the same performance as the time-independent implementation of ContinuousWave; in fact, this is how ContinuousWave is implemented for time dependent arguments. However, this manual hack is discouraged, as ContinuousWave provides a simpler interface, that clearly expresses its intend and is potentially faster for constant driving amplitudes. For normalisable wave functions both Coherent and ContinuousWave are practically interchangeable, but it is still recommended to use Coherent, as in the last example, to express clear intend.","category":"page"},{"location":"examples/CoherentInput/#Accessing-the-input-cavity","page":"Displaced States","title":"Accessing the input cavity","text":"","category":"section"},{"location":"examples/CoherentInput/","page":"Displaced States","title":"Displaced States","text":"All previous examples eliminate the input cavity in favour for better performance. If one wants to access the precise population of the input mode however, we can save the coherent number distribution in ArbitraryState. This needs a little boilerplating","category":"page"},{"location":"examples/CoherentInput/","page":"Displaced States","title":"Displaced States","text":"N_cutoff = ceil(Int, tf * α^2 + 5 * √(tf) * α)\namplitudes = coherentstate(FockBasis(N_cutoff), √(tf) * α).data\nstate = ArbitraryState(amplitudes)","category":"page"},{"location":"examples/CoherentInput/","page":"Displaced States","title":"Displaced States","text":"First, a cutoff is choosen, so that the coherent state fits reasonably well in the Hilbert space (here sum(abs2, amplitudes) == 0.9999...). Then we use QuantumOptics.jl coherentstate to get the amplitudes and create our input state. Notice, we choose an amplitude √(tf) * α and not α, since we use HardBoxMode as our mode function, which introduces a normalisation factor 1/√(tf). Here it is important that we use a proper normalisable mode, since we need to access gᵢ. The simulation now becomes","category":"page"},{"location":"examples/CoherentInput/","page":"Displaced States","title":"Displaced States","text":"problem = WaveguideProblem(TwoLevelChain(1), WavePacket(HardBoxMode(σ = tf), state), tf)\nts, ρs_full = @btime solve($problem) # 1.300 s (282065 allocations: 560.71 MiB)","category":"page"},{"location":"examples/CoherentInput/","page":"Displaced States","title":"Displaced States","text":"and is the significant slowest version. Yet, we now may observe the input cavity and its entanglement with the quantum system or output cavity","category":"page"},{"location":"examples/CoherentInput/","page":"Displaced States","title":"Displaced States","text":"basis(ρs_full[end]) # [Fock(cutoff=17) ⊗ NLevel(N=2)]","category":"page"},{"location":"examples/CoherentInput/","page":"Displaced States","title":"Displaced States","text":"Lastly, we verify that this approach yields the correct result for the atom's state","category":"page"},{"location":"examples/CoherentInput/","page":"Displaced States","title":"Displaced States","text":"ρs_atoms = [ptrace(ρ, 1) for ρ ∈ ρs_full]\nall(tracedistance.(ρs_atoms, ρs) .< 2e-5) # true","category":"page"},{"location":"examples/CustomSystems/#Defining-the-quantum-system","page":"Defining the quantum system","title":"Defining the quantum system","text":"","category":"section"},{"location":"examples/CustomSystems/","page":"Defining the quantum system","title":"Defining the quantum system","text":"Each quantum system used in a simulation consists of three ingredients:","category":"page"},{"location":"examples/CustomSystems/","page":"Defining the quantum system","title":"Defining the quantum system","text":"The generator of the system's dynamic in the absence of the input cavity.\nThe operator that couples the external bosonic modes to the quantum system.\nThe initial state of the system.","category":"page"},{"location":"examples/CustomSystems/","page":"Defining the quantum system","title":"Defining the quantum system","text":"These parameters are simply passed to WaveguideProblem as an array","category":"page"},{"location":"examples/CustomSystems/","page":"Defining the quantum system","title":"Defining the quantum system","text":"WaveguideProblem([hamiltonian, σ, ψ₀], WavePacket(...), t)","category":"page"},{"location":"examples/CustomSystems/","page":"Defining the quantum system","title":"Defining the quantum system","text":"or, for dissipative dynamic,","category":"page"},{"location":"examples/CustomSystems/","page":"Defining the quantum system","title":"Defining the quantum system","text":"WaveguideProblem([hamiltonian, [Ls], σ, ψ₀], WavePacket(...), t)","category":"page"},{"location":"examples/CustomSystems/","page":"Defining the quantum system","title":"Defining the quantum system","text":"It should be noted that, right now, the initial state of the quantum system cannot be entangled with the input or output modes.","category":"page"},{"location":"examples/CustomSystems/","page":"Defining the quantum system","title":"Defining the quantum system","text":"As a practical example we consider a two level atom, which we model as a spin-1/2 particle. The atom has a dephasing δ, spontaneous decay rate Γ, starts in the excited state, and is coupled to the bosonic modes via the interaction Hamiltonian","category":"page"},{"location":"examples/CustomSystems/","page":"Defining the quantum system","title":"Defining the quantum system","text":"H = - fraci2 bigσ^+ b(0) - σ^- b^dagger(0)","category":"page"},{"location":"examples/CustomSystems/","page":"Defining the quantum system","title":"Defining the quantum system","text":"The operator that couples the quantum system to the bosons is defined to be the operator next to the b^dagger(0) term, i.e., σ^- here. To simulate this system we use QuantumOptics.jl's spin operators","category":"page"},{"location":"examples/CustomSystems/","page":"Defining the quantum system","title":"Defining the quantum system","text":"basis  = SpinBasis(1//2)\nσ⁺, σ⁻ = sigmap(basis), sigmam(basis)\n\nδ, Γ = rand(2)\n\nH  = δ * σ⁺ * σ⁻\nL  = √(Γ) * σ⁻\nψ₀ = spinup(basis)","category":"page"},{"location":"examples/CustomSystems/","page":"Defining the quantum system","title":"Defining the quantum system","text":"Then, we simply create the WaveguideProblem","category":"page"},{"location":"examples/CustomSystems/","page":"Defining the quantum system","title":"Defining the quantum system","text":"WaveguideProblem([H, [L], σ⁻, ψ₀], WavePacket(...), t)","category":"page"},{"location":"examples/CustomSystems/","page":"Defining the quantum system","title":"Defining the quantum system","text":"By default ChiralWaveguide.jl assumes that the interaction Hamiltonian between bosons and the quantum system has coupling rate 1. To change this you should rescale all other rates (including the bosons' rates) by the respective coupling rate you wish for. Lastly, nothing is special about the spin-1/2 representation we choose here. In fact the exported systems TwoLevelChain and DissipativeLambdaChain use NLevelBasis instead and every basis in the QuantumOptics.jl ecosystem works.","category":"page"},{"location":"examples/SinglePhotonScattering/#Scattering-of-a-single-photon","page":"Scattering of a single photon","title":"Scattering of a single photon","text":"","category":"section"},{"location":"examples/SinglePhotonScattering/","page":"Scattering of a single photon","title":"Scattering of a single photon","text":"In this example, we analyse the scattering of a single photon in a Gaussian mode at a single chiral two level atom. For comparison, we reconstruct the second figure from the paper from A. Kiilerich, K. Mølmer, on which this package is based. For this problem, the scatter problem can be solved analytically and the the outgoing mode function is","category":"page"},{"location":"examples/SinglePhotonScattering/","page":"Scattering of a single photon","title":"Scattering of a single photon","text":"\tpsi_mathrmout(t) propto expleft(-frac(t - τ)^22right)\n\t\t- sqrtfracπ2 expleft(- fract - τ2 + frac18 right)\n\t\t  mathrmerfcleft(frac-2(t - τ) + 12 sqrt2right)","category":"page"},{"location":"examples/SinglePhotonScattering/","page":"Scattering of a single photon","title":"Scattering of a single photon","text":"where the incoming Gauss mode has a width of sigma = 1 and it's center of mass arrives at the atom at time tau.","category":"page"},{"location":"examples/SinglePhotonScattering/","page":"Scattering of a single photon","title":"Scattering of a single photon","text":"We start the simulation by loading the necessary packages","category":"page"},{"location":"examples/SinglePhotonScattering/","page":"Scattering of a single photon","title":"Scattering of a single photon","text":"using ChiralWaveguide, Plots, SpecialFunctions # SpecialFunctions has erfc","category":"page"},{"location":"examples/SinglePhotonScattering/","page":"Scattering of a single photon","title":"Scattering of a single photon","text":"defining the incoming single photon wave packet","category":"page"},{"location":"examples/SinglePhotonScattering/","page":"Scattering of a single photon","title":"Scattering of a single photon","text":"wavepacket = WavePacket(GaussMode(τ = 4.0), Fock(1))","category":"page"},{"location":"examples/SinglePhotonScattering/","page":"Scattering of a single photon","title":"Scattering of a single photon","text":"and create the final mode we want to observe (notice the mode does not have to be normalised)","category":"page"},{"location":"examples/SinglePhotonScattering/","page":"Scattering of a single photon","title":"Scattering of a single photon","text":"function outputModeFunction(t, τ = 4.0)\n\tabs(t-τ) > 10 && return 0.0 # cutoff for numerical stability\n\treturn exp(-(t - τ)^2/2) -\n\t\t√(π/2) * exp(-(t - τ) / 2 + 1 / 8) * erfc((-2(t - τ) + 1) / (2 * √(2)))\nend\n\noutputMode = Mode(t -> outputModeFunctionFig2(t, 1.0, 4.0))","category":"page"},{"location":"examples/SinglePhotonScattering/","page":"Scattering of a single photon","title":"Scattering of a single photon","text":"Notice that we did not provide the coupling rates for the outputMode ourself but let the constructor of Mode take care of that, which numerically solves the integrals for us.","category":"page"},{"location":"examples/SinglePhotonScattering/","page":"Scattering of a single photon","title":"Scattering of a single photon","text":"We now have every constituents to define","category":"page"},{"location":"examples/SinglePhotonScattering/","page":"Scattering of a single photon","title":"Scattering of a single photon","text":"problem = WaveguideProblem(TwoLevelChain(1), wavepacket, outputMode, 13.0)","category":"page"},{"location":"examples/SinglePhotonScattering/","page":"Scattering of a single photon","title":"Scattering of a single photon","text":"and solve the scattering problem","category":"page"},{"location":"examples/SinglePhotonScattering/","page":"Scattering of a single photon","title":"Scattering of a single photon","text":"ts, ρs = solve(problem)","category":"page"},{"location":"examples/SinglePhotonScattering/","page":"Scattering of a single photon","title":"Scattering of a single photon","text":"This already completes the simulation and we can now define the observables of interest. The density matrices ρs at time ts are ordinary QuantumOptics.jl operators and we can therefore utilise every tool within QuantumOptics.jl.","category":"page"},{"location":"examples/SinglePhotonScattering/","page":"Scattering of a single photon","title":"Scattering of a single photon","text":"basis  = ρs[end].basis_l\nn̂ᵢ, n̂ₒ = number(basis.bases[1]), number(basis.bases[3])\nσ⁺σ⁻   = transition(NLevelBasis(2), 2, 2)","category":"page"},{"location":"examples/SinglePhotonScattering/","page":"Scattering of a single photon","title":"Scattering of a single photon","text":"Finally, we draw the populations of the input cavity, the excited state of the atom, and the output cavity, resulting in the same figure as from the paper:","category":"page"},{"location":"examples/SinglePhotonScattering/","page":"Scattering of a single photon","title":"Scattering of a single photon","text":"plot( ts, expect(1, n̂ᵢ,   ρs) .|> real, label = \"⟨n̂_i⟩\",  line = (2, :red))\nplot!(ts, expect(2, σ⁺σ⁻, ρs) .|> real, label = \"⟨σ⁺σ⁻⟩\", line = (2, :green, :dot))\nplot!(ts, expect(3, n̂ₒ,   ρs) .|> real, label = \"⟨n̂_o⟩\",  line = (2, :blue,  :dash))\nplot!(size = (600, 200), legend = :right)","category":"page"},{"location":"examples/SinglePhotonScattering/","page":"Scattering of a single photon","title":"Scattering of a single photon","text":"(Image: Populations)","category":"page"},{"location":"examples/SinglePhotonScattering/","page":"Scattering of a single photon","title":"Scattering of a single photon","text":"","category":"page"},{"location":"examples/SinglePhotonScattering/","page":"Scattering of a single photon","title":"Scattering of a single photon","text":"Entire Script","category":"page"},{"location":"examples/SinglePhotonScattering/","page":"Scattering of a single photon","title":"Scattering of a single photon","text":"using ChiralWaveguide, Plots, SpecialFunctions # SpecialFunctions has erfc\n\nwavepacket = WavePacket(GaussMode(τ = 4.0), Fock(1))\n\nfunction outputModeFunction(t, τ)\n\tabs(t-τ) > 10 && return 0.0 # cutoff for numerical stability\n\treturn exp(-(t - τ)^2/2) -\n\t\t√(π/2) * exp(-(t - τ) / 2 + 1 / 8) * erfc((-2(t - τ) + 1) / (2 * √(2)))\nend\n\noutputMode = Mode(t -> outputModeFunction(t, 4.0))\n\nproblem = WaveguideProblem(TwoLevelChain(1), wavepacket, outputMode, 13.0)\nts, ρs = solve(problem)\n\nbasis  = ρs[end].basis_l\nn̂ᵢ, n̂ₒ = number(basis.bases[1]), number(basis.bases[3])\nσ⁺σ⁻   = transition(NLevelBasis(2), 2, 2)\n\nplot( ts, expect(1, n̂ᵢ,   ρs) .|> real, label = \"⟨n̂_i⟩\",  line = (2, :red))\nplot!(ts, expect(2, σ⁺σ⁻, ρs) .|> real, label = \"⟨σ⁺σ⁻⟩\", line = (2, :green, :dot))\nplot!(ts, expect(3, n̂ₒ,   ρs) .|> real, label = \"⟨n̂_o⟩\",  line = (2, :blue,  :dash))\nplot!(size = (600, 200), legend = :right)","category":"page"},{"location":"#Simulating-1d-chiral-quantum-systems","page":"Simulating 1d chiral quantum systems","title":"Simulating 1d chiral quantum systems","text":"","category":"section"},{"location":"","page":"Simulating 1d chiral quantum systems","title":"Simulating 1d chiral quantum systems","text":"ChiralWaveguide.jl is a simulation package for arbitrary chiral quantum systems, coupled to a continuous bosonic background with linear dispersion, for example photons. This package aims to deliver three promises:","category":"page"},{"location":"","page":"Simulating 1d chiral quantum systems","title":"Simulating 1d chiral quantum systems","text":"A lightweight, descriptive interface. Scattering a Gaussian Fock state on a chiral atom is no more complicated than\nproblem = WaveguideProblem(TwoLevelChain(1), WavePacket(GaussMode(), Fock(1)), 10.0)\nts, ρs = solve(problem)\nPerformance. ChiralWaveguide.jl performs suitable transformations on the Hilbert space, like the Mollow transform, to reduce the effective Hilbert space size as much as possible.\nInteroperability and a native Julia experience. ChiralWaveguide is build entirely around the QuantumOptics.jl package. In fact, getting the population of the excited state in the previous example is no more complicated than\nσ⁺σ⁻ = transition(NLevelBasis(2), 2, 2)\nexpect(2, σ⁺σ⁻, ρs)\nAdditionally, the syntax is heavily inspired by the SciML ecosystem.","category":"page"}]
}
