# Objective: the composite normalized error (CNE) of Bellamy et al. (2026),
#     CNE_β = w_mlr / max|MLR_exp| · MAE(MLR) + w_m / max(m_exp) · MAE(m)
#     CNE   = Σ_β w_β CNE_β          (w_β normalized to sum to 1)
# over the compared points inside fit_range. 0 is a perfect fit.

const PENALTY = 1e10

struct FitProblem{N, M}
	scheme::Scheme
	datasets::Vector{Dataset}
	weights::Vector{Float64}
	w_mass::Float64
	w_mlr::Float64
	opts::SolverOptions
	tsave::Vector{Vector{Float64}}      # compared times up to the last one inside fit_range
	fidx::Vector{Vector{Int}}           # indices of those inside fit_range
end

function FitProblem(scheme::Scheme, datasets::Vector{Dataset}, obj::AbstractDict, opts::SolverOptions)
	all(has_data, datasets) || cfgerror("fitting needs a data file for every [[dataset]]")
	w = [d.weight for d in datasets]
	sum(w) > 0 || cfgerror("dataset weights add up to zero")
	w ./= sum(w)
	wm = as_float(obj["mass_weight"], "mass_weight")
	wr = as_float(obj["mlr_weight"], "mlr_weight")
	(wm >= 0 && wr >= 0 && wm + wr > 0) || cfgerror("objective weights must be >= 0 and not both zero")
	tot = wm + wr
	N, M = length(scheme.components), length(scheme.reactions)
	fidx = [findall(d.infit) for d in datasets]
	tsave = [d.t[1:last(i)] for (d, i) in zip(datasets, fidx)]
	return FitProblem{N, M}(scheme, datasets, w, wm / tot, wr / tot, opts, tsave, fidx)
end

kinetics(fp::FitProblem{N, M}, x) where {N, M} = build_kinetics(fp.scheme, x, Val(N), Val(M))

"""
    objective(fp, x) -> Float64

CNE for the optimizer vector `x` (lower is better). Infeasible yields or
initial fractions, and failed solves, return a large penalty.
"""
function objective(fp::FitProblem, x::AbstractVector{Float64})
	kin = kinetics(fp, x)
	kin === nothing && return PENALTY
	total = 0.0
	for (k, d) in enumerate(fp.datasets)
		ts = fp.tsave[k]
		u, ok = simulate_run(kin, d.prog, ts, fp.opts)
		ok || return PENALTY
		total += fp.weights[k] * cne(d, u, ts, fp.fidx[k], fp.w_mass, fp.w_mlr)
	end
	return isfinite(total) ? total : PENALTY
end

function cne(d::Dataset, u, ts, idx, w_mass, w_mlr)
	mass = [sum(v) for v in u]
	mlr = mass_loss_rate!(similar(mass), mass, ts)
	return cne_from_curves(d, mass, mlr, idx, w_mass, w_mlr)
end

function cne_from_curves(d::Dataset, mass, mlr, idx, w_mass, w_mlr)
	sm = 0.0
	sr = 0.0
	@inbounds for i in idx
		sm += abs(mass[i] - d.mass[i])
		sr += abs(mlr[i] - d.mlr[i])
	end
	n = length(idx)
	return w_mlr / d.mlr_max * (sr / n) + w_mass / d.mass_max * (sm / n)
end

"""
    model_curves(kin, d; opts=TIGHT)

Model mass, mass-loss rate and component masses at every compared point of
dataset `d`, solved with tight tolerances.
"""
function model_curves(kin::Kinetics{N}, d::Dataset; opts::SolverOptions = TIGHT) where {N}
	u, ok = simulate_run(kin, d.prog, d.t, opts)
	ok || error("the model could not be solved for dataset $(d.name)")
	mass = [sum(v) for v in u]
	mlr = mass_loss_rate!(similar(mass), mass, d.t)
	comps = [[v[i] for v in u] for i in 1:N]
	return mass, mlr, comps
end

function cne_tight(kin::Kinetics, d::Dataset, w_mass, w_mlr)
	mass, mlr, _ = model_curves(kin, d)
	return cne_from_curves(d, mass, mlr, findall(d.infit), w_mass, w_mlr)
end
