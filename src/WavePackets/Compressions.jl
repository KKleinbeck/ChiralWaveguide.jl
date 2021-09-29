struct Compression
	comp::Function
	decomp::Function
	inv_jacobian::Function

	σ::Float64
	μ::Float64
	function Compression(c, d, i; σ = 1.0, μ = 0.0)
		σ ≤ 0.0 && error("σ needs to be positive")
		new(c, d, i, σ, μ)
	end
end

Algebraic(; σ = 1.0, μ = 0.0) = Compression(
	t ->   t / (1.0 + abs(t)),
	τ ->   τ / (1.0 - abs(τ)),
	τ -> 1.0 / (1.0 - abs(τ))^2, # == (1.0 + abs( τ / (1.0 - abs(τ)) ) )^2
	σ = σ, μ = μ
)

Exponential(; σ = 1.0, μ = 0.0) = Compression(
	t -> tanh(t),
	τ -> atanh(τ),
	τ -> 1.0 / (1.0 - τ^2), # == cosh(atanh(τ))^2
	σ = σ, μ = μ
)

Trigonometric(; σ = 1.0, μ = 0.0) = Compression(
	t -> 2.0 * atan(t) / π,
	τ -> tan(π*τ/2.0),
	τ -> π * (1.0 + tan(π*τ/2.0)^2) / 2.0,
	σ = σ, μ = μ
)

Compressions = Dict(
	:algebraic     => Algebraic(),
	:exponential   => Exponential(),
	:trigonometric => Trigonometric()
)
