# Mass-loss rate from measured mass: a cubic P-spline (equally spaced knots,
# second-order difference penalty; Eilers and Marx) fitted to mass against
# time, with the smoothing weight λ picked by the V-curve (Frasso and Eilers).
# The derivative of the fitted spline is the mass-loss rate.

struct PSplineFit
	xl::Float64
	h::Float64      # knot spacing
	nseg::Int
	coef::Vector{Float64}
	lambda::Float64
end

# Cubic B-spline values and first derivatives (per unit u) on one segment.
@inline function _cubic_basis(u)
	v = 1 - u
	b = (v^3 / 6, (3u^3 - 6u^2 + 4) / 6, (-3u^3 + 3u^2 + 3u + 1) / 6, u^3 / 6)
	d = (-v^2 / 2, (3u^2 - 4u) / 2, (-3u^2 + 2u + 1) / 2, u^2 / 2)
	return b, d
end

@inline function _segment(x, xl, h, nseg)
	s = (x - xl) / h
	k = clamp(floor(Int, s), 0, nseg - 1)
	return k, s - k
end

function _basis_matrix(x, xl, h, nseg)
	n = length(x)
	I = Vector{Int}(undef, 4n)
	J = Vector{Int}(undef, 4n)
	V = Vector{Float64}(undef, 4n)
	@inbounds for i in 1:n
		k, u = _segment(x[i], xl, h, nseg)
		b, _ = _cubic_basis(u)
		for m in 1:4
			idx = 4(i - 1) + m
			I[idx] = i
			J[idx] = k + m
			V[idx] = b[m]
		end
	end
	return sparse(I, J, V, n, nseg + 3)
end

function _diff_matrix(nb, order)
	D = sparse(1.0I, nb, nb)
	for _ in 1:order
		D = D[2:end, :] - D[1:(end-1), :]
	end
	return D
end

"""
    pspline_fit(x, y; nseg=500, lambda=:auto, lambda_bounds=(0.1, 1e4), n_lambda=200)

Cubic P-spline fit of `y` against `x`. With `lambda = :auto` the smoothing
weight is the V-curve choice over a log grid within `lambda_bounds`.
"""
function pspline_fit(x::AbstractVector, y::AbstractVector; nseg::Int = 500, lambda = :auto,
	lambda_bounds = (0.1, 1e4), n_lambda::Int = 200)
	xl, xr = extrema(x)
	h = (xr - xl) / nseg
	B = _basis_matrix(x, xl, h, nseg)
	D = _diff_matrix(nseg + 3, 2)
	BtB = B' * B
	Bty = B' * y
	DtD = D' * D
	solve_for(λ) = cholesky(Symmetric(BtB + λ * DtD)) \ Bty

	if lambda === :auto
		grid = exp10.(range(log10(lambda_bounds[1]), log10(lambda_bounds[2]); length = n_lambda))
		lr = similar(grid)
		lp = similar(grid)
		for (k, λ) in enumerate(grid)
			a = solve_for(λ)
			lr[k] = log10(sum(abs2, y - B * a))
			lp[k] = log10(sum(abs2, D * a))
		end
		ok = isfinite.(lr) .& isfinite.(lp)
		count(ok) >= 2 || error("P-spline V-curve: too few finite points")
		g, lr, lp = grid[ok], lr[ok], lp[ok]
		dist = hypot.(diff(lr), diff(lp))
		k = argmin(dist)
		lambda = sqrt(g[k] * g[k+1])
	end
	λ = Float64(lambda)
	return PSplineFit(xl, h, nseg, solve_for(λ), λ)
end

@inline function spline_value(p::PSplineFit, x::Real)
	k, u = _segment(x, p.xl, p.h, p.nseg)
	b, _ = _cubic_basis(u)
	@inbounds return b[1] * p.coef[k+1] + b[2] * p.coef[k+2] + b[3] * p.coef[k+3] + b[4] * p.coef[k+4]
end

@inline function spline_slope(p::PSplineFit, x::Real)
	k, u = _segment(x, p.xl, p.h, p.nseg)
	_, d = _cubic_basis(u)
	@inbounds return (d[1] * p.coef[k+1] + d[2] * p.coef[k+2] + d[3] * p.coef[k+3] + d[4] * p.coef[k+4]) / p.h
end

function pspline_eval(p::PSplineFit, x::AbstractVector; derivative::Bool = false)
	out = similar(x, Float64)
	@inbounds for i in eachindex(x)
		k, u = _segment(x[i], p.xl, p.h, p.nseg)
		b, d = _cubic_basis(u)
		w = derivative ? d : b
		acc = 0.0
		for m in 1:4
			acc += w[m] * p.coef[k+m]
		end
		out[i] = derivative ? acc / p.h : acc
	end
	return out
end
