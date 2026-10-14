# The kinetic model. State u = mass of each component per unit reference mass.
# Reaction j consumes its reactant a at
#     r_j = A_j exp(-E_j / (R T)) u_a^n_j
# and makes ν_jp r_j of each solid product p; the rest, (1 - Σ_p ν_jp) r_j,
# leaves as gas. Sample mass = Σ u. This is the FDS solid-phase reaction form
# with the mass normalized by the reference mass. The right-hand side and its
# Jacobian are written out by hand on static arrays, so a solve does not
# allocate.
#
# The model mass-loss rate is the time derivative of the solved sample mass
# (three-point differences), as in the original scripts. Summing the gas
# release Σ_j (1 - Σ_p ν_jp) r_j at the saved states is equal in exact
# arithmetic but multiplies the solver's error in a fast reaction's reactant
# by its rate constant.

struct Kinetics{N, M, L}
	A::SVector{M, Float64}
	EoR::SVector{M, Float64}       # E / R, K
	order::SVector{M, Float64}
	reactant::SVector{M, Int}
	S::SMatrix{N, M, Float64, L}   # net change of component i per unit r_j
	gas::SVector{M, Float64}       # gas fraction of reaction j
	y0::SVector{N, Float64}
end

struct ODEParams{K <: Kinetics}
	kin::K
	prog::TempProgram
end

# First order: k·y, linear and smooth through zero, so a reactant that a stiff
# solve leaves slightly negative is pulled back to zero. Other orders use
# max(y, 0)^n, which needs y ≥ 0.
@inline _conc(y, n) = n == 1.0 ? y : (y > 0 ? y^n : 0.0)
@inline _dconc(y, n) = n == 1.0 ? 1.0 : (y > 0 ? n * y^(n - 1) : 0.0)

@inline function rates(u, kin::Kinetics{N, M}, T) where {N, M}
	invT = 1 / T
	return SVector{M, Float64}(ntuple(Val(M)) do j
		@inbounds kin.A[j] * exp(-kin.EoR[j] * invT) * _conc(u[kin.reactant[j]], kin.order[j])
	end)
end

function rhs(u, p::ODEParams, t)
	T = temperature(p.prog, t)
	return p.kin.S * rates(u, p.kin, T)
end

function jac(u, p::ODEParams{<:Kinetics{N, M}}, t) where {N, M}
	kin = p.kin
	invT = 1 / temperature(p.prog, t)
	dr = SVector{M, Float64}(ntuple(Val(M)) do j
		@inbounds kin.A[j] * exp(-kin.EoR[j] * invT) * _dconc(u[kin.reactant[j]], kin.order[j])
	end)
	# D[j, a] = ∂r_j/∂u_a, nonzero only for a = reactant of j
	D = SMatrix{M, N, Float64}(ntuple(Val(M * N)) do idx
		j = (idx - 1) % M + 1
		a = (idx - 1) ÷ M + 1
		@inbounds kin.reactant[j] == a ? dr[j] : 0.0
	end)
	return kin.S * D
end

function tgrad(u, p::ODEParams, t)
	T = temperature(p.prog, t)
	dTdt = temperature_slope(p.prog, t)
	r = rates(u, p.kin, T)
	return p.kin.S * (r .* p.kin.EoR .* (dTdt / T^2))
end

const ODE_FUNCTION = ODEFunction{false}(rhs; jac = jac, tgrad = tgrad)

"""
    mass_loss_rate!(mlr, m, t)

-dm/dt at every point from second-order three-point differences (one-sided
at the ends), for any spacing of `t`.
"""
function mass_loss_rate!(mlr, m, t)
	n = length(m)
	@inbounds for i in 1:n
		j = clamp(i, 2, n - 1)
		h1 = t[j] - t[j-1]
		h2 = t[j+1] - t[j]
		a, b, c = m[j-1], m[j], m[j+1]
		if i == j
			d = (h1^2 * c - h2^2 * a + (h2^2 - h1^2) * b) / (h1 * h2 * (h1 + h2))
		elseif i < j
			d = (-(2h1 + h2) * h2 * a + (h1 + h2)^2 * b - h1^2 * c) / (h1 * h2 * (h1 + h2))
		else
			d = (h2^2 * a - (h1 + h2)^2 * b + (h1 + 2h2) * h1 * c) / (h1 * h2 * (h1 + h2))
		end
		mlr[i] = -d
	end
	return mlr
end

# Values of the scheme's scalars for one optimizer vector x. Returns `nothing`
# when the yields of a reaction or the initial fractions are not feasible.
function build_kinetics(s::Scheme, x, ::Val{N}, ::Val{M}) where {N, M}
	A = MVector{M, Float64}(undef)
	EoR = MVector{M, Float64}(undef)
	order = MVector{M, Float64}(undef)
	reactant = MVector{M, Int}(undef)
	S = zeros(MMatrix{N, M, Float64})
	gas = MVector{M, Float64}(undef)
	@inbounds for (j, r) in enumerate(s.reactions)
		A[j] = exp10(slotvalue(r.log10A, x))
		EoR[j] = slotvalue(r.E, x) * 1e3 / R_GAS
		order[j] = slotvalue(r.order, x)
		reactant[j] = r.reactant
		S[r.reactant, j] = -1.0
		total = 0.0
		for (p, y) in zip(r.products, r.yields)
			v = slotvalue(y, x)
			S[p, j] += v
			total += v
		end
		total > 1.0 + 1e-12 && return nothing
		gas[j] = max(1.0 - total, 0.0)
	end
	y0 = MVector{N, Float64}(undef)
	acc = 0.0
	@inbounds for i in 1:N
		y0[i] = i == s.rest ? 0.0 : slotvalue(s.initial[i], x)
		acc += y0[i]
	end
	if s.rest != 0
		acc > 1.0 + 1e-12 && return nothing
		y0[s.rest] = max(1.0 - acc, 0.0)
	end
	return Kinetics{N, M, N * M}(SVector(A), SVector(EoR), SVector(order), SVector(reactant),
		SMatrix(S), SVector(gas), SVector(y0))
end

# Out-of-domain guard (a component going clearly negative), as in the
# original scripts.
@inline function _isoutofdomain(u, p, t)
	@inbounds for v in u
		v < -1e-5 && return true
	end
	return false
end

const PRIMARY_SOLVER = Rodas5P(autodiff = AutoFiniteDiff(), max_jac_age = 1)
const FALLBACK_SOLVER = Rosenbrock23(autodiff = AutoFiniteDiff())

struct SolverOptions
	abstol::Float64
	reltol::Float64
	maxiters::Int
end

SolverOptions(d::AbstractDict) = SolverOptions(as_float(d["abstol"], "abstol"),
	as_float(d["reltol"], "reltol"), as_int(d["maxiters"], "maxiters"))

const TIGHT = SolverOptions(1e-9, 1e-7, 1_000_000)
const QUIET = DEVerbosity(SciMLLogging.None())

function _init(prob, alg, ts, o::SolverOptions)
	return init(prob, alg; saveat = ts, abstol = o.abstol, reltol = o.reltol,
		maxiters = o.maxiters, save_everystep = false, dense = false,
		isoutofdomain = _isoutofdomain, verbose = QUIET)
end

"""
    simulate_run(kin, prog, ts, opts) -> (u, ok)

Solve the scheme on one temperature program and return the state at the
times `ts` (which start at or after `prog.t0`). If the main solver fails, a
lower-order Rosenbrock method is tried before giving up.
"""
function simulate_run(kin::Kinetics, prog::TempProgram, ts::Vector{Float64}, opts::SolverOptions)
	prob = ODEProblem{false}(ODE_FUNCTION, kin.y0, (prog.t0, ts[end]), ODEParams(kin, prog))
	for alg in (PRIMARY_SOLVER, FALLBACK_SOLVER)
		u = _try_solve(prob, alg, ts, opts)
		u === nothing || return u, true
	end
	return nothing, false
end

function _try_solve(prob, alg, ts, opts)
	try
		integ = _init(prob, alg, ts, opts)
		solve!(integ)
		ok = successful_retcode(integ.sol.retcode) && length(integ.sol.u) == length(ts)
		return ok ? integ.sol.u : nothing
	catch e
		e isa InterruptException && rethrow()
		return nothing
	end
end
