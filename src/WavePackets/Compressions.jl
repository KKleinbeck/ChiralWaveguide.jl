struct Compression
	comp::Function
	decomp::Function
	inv_jacobian::Function
	scale::Float64
end

Algebraic(s = 1.0) = Compression(
	t ->   t / (1.0 + abs(t)),
	τ ->   τ / (1.0 - abs(τ)),
	τ -> 1.0 / (1.0 - abs(τ))^2, # == (1.0 + abs( τ / (1.0 - abs(τ)) ) )^2
	s
)

Exponential(s = 1.0) = Compression(
	t -> tanh(t),
	τ -> atanh(τ),
	τ -> 1.0 / (1.0 - τ^2), # == cosh(atanh(τ))^2
	s
)

Trigonometric(s = 1.0) = Compression(
	t -> 2.0 * atan(t) / π,
	τ -> tan(π*τ/2.0),
	τ -> π * (1.0 + tan(π*τ/2.0)^2) / 2.0,
	s
)

Compressions = Dict(
	:algebraic     => Algebraic(),
	:exponential   => Exponential(),
	:trigonometric => Trigonometric()
)
