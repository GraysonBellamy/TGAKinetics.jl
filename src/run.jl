# The three commands: check, fit, simulate.

function prepare(path)
	cfg = load_config(path)
	scheme = parse_scheme(cfg)
	raw = get(cfg, "dataset", nothing)
	(raw isa AbstractVector && !isempty(raw)) || cfgerror("no [[dataset]] entries")
	datasets = [load_dataset(cfg, e, k) for (k, e) in enumerate(raw)]
	allunique(d.name for d in datasets) || cfgerror("dataset names must be unique")
	return cfg, scheme, datasets
end

function describe_datasets(io, datasets)
	for d in datasets
		println(io, "Dataset ", d.name, has_data(d) ? " ($(basename(d.file)))" : " (set heating program)")
		if has_data(d)
			ipk = argmax(d.mlr)
			@printf(io, "  reference mass %.6g at %.1f K; %d compared points, %d inside fit_range (%.1f-%.1f K)\n",
				d.ref_mass, d.ref_temperature, length(d.t), count(d.infit),
				minimum(d.T[d.infit]), maximum(d.T[d.infit]))
			@printf(io, "  peak mass-loss rate %.4e 1/s at %.1f K; final mass %.4f\n",
				d.mlr[ipk], d.T[ipk], d.mass[end])
			isnan(d.lambda) || @printf(io, "  P-spline smoothing weight λ = %.4g\n", d.lambda)
		else
			@printf(io, "  %.1f-%.1f K\n", d.T[1], d.T[end])
		end
		@printf(io, "  weight %.4g\n", d.weight)
	end
end

"""
    check(path)

Read the run file and data, report what was found, write the processed data
(`data_<name>.csv`) to the output folder, and time the objective to estimate
how long `fit` will take. Run this before a long fit.
"""
function check(path::AbstractString; io::IO = stdout)
	cfg, scheme, datasets = prepare(path)
	outdir = output_dir(cfg)
	describe(io, scheme)
	describe_datasets(io, datasets)
	for d in datasets
		has_data(d) && write_processed(joinpath(outdir, "data_$(safe_name(d.name)).csv"), d)
	end
	println(io, "Processed data written to ", outdir)
	all(has_data, datasets) || return nothing

	fp = FitProblem(scheme, datasets, objective_settings(cfg), SolverOptions(solver_settings(cfg)))
	n = n_free(scheme)
	if n == 0
		@printf(io, "No free parameters. CNE = %.6e\n", objective(fp, Float64[]))
		return nothing
	end
	mid = (scheme.lower .+ scheme.upper) ./ 2
	objective(fp, mid)  # compile
	rng = Xoshiro(0)
	samples = [scheme.lower .+ rand(rng, n) .* (scheme.upper .- scheme.lower) for _ in 1:20]
	times = [(@elapsed objective(fp, x)) for x in samples]
	per_eval = sort(times)[10]
	opt = DEOptions(optimizer_settings(cfg))
	evals = opt.runs * opt.population * (opt.generations + 1)
	@printf(io, "One evaluation (all datasets) takes about %.1f ms on one thread.\n", 1e3 * per_eval)
	@printf(io, "Population %d, %d generations, %d run(s): up to %d evaluations, roughly %s with %d threads\n",
		opt.population, opt.generations, opt.runs, evals,
		format_duration(per_eval * evals / Threads.nthreads()), Threads.nthreads())
	println(io, "(timed at random points inside the bounds; near a good fit an evaluation can cost more)")
	return nothing
end

"""
    fit(path) -> NamedTuple

Fit the free parameters of the run file. Writes to the output folder:
`parameters.csv`, `fitted.toml` (a run file with the fitted values, ready for
`simulate`), `fit_<name>.csv` (measured and model curves), `convergence.csv`
and `summary.txt`. While running, the best parameters so far are kept in
`best_so_far.csv`.
"""
function fit(path::AbstractString; io::IO = stdout)
	cfg, scheme, datasets = prepare(path)
	outdir = output_dir(cfg)
	n_free(scheme) > 0 || cfgerror("nothing to fit: every parameter has a fixed value (use simulate)")
	fp = FitProblem(scheme, datasets, objective_settings(cfg), SolverOptions(solver_settings(cfg)))
	opt = DEOptions(optimizer_settings(cfg))
	describe(io, scheme)
	describe_datasets(io, datasets)
	@printf(io, "Fitting %d parameters: population %d, up to %d generations, %d threads\n",
		n_free(scheme), opt.population, opt.generations, Threads.nthreads())
	flush(io)

	checkpoint = joinpath(outdir, "best_so_far.csv")
	results = []
	best_f = Inf
	for run in 1:opt.runs
		opt.runs > 1 && println(io, "\nRun $run of $(opt.runs) (seed $(opt.seed + run - 1))")
		r = de1220(x -> objective(fp, x), scheme.lower, scheme.upper, opt; run = run, io = io,
			on_report = (g, x, f) -> f <= best_f && write_parameters(checkpoint, scheme, x))
		println(io, "Stopped: ", r.stop)
		push!(results, r)
		best_f = min(best_f, r.f)
	end
	res = results[argmin([r.f for r in results])]
	opt.runs > 1 && write_runs(joinpath(outdir, "runs.csv"), results)

	kin = kinetics(fp, res.x)
	per = [cne_tight(kin, d, fp.w_mass, fp.w_mlr) for d in datasets]
	total = sum(fp.weights .* per)

	write_parameters(joinpath(outdir, "parameters.csv"), scheme, res.x)
	write_fitted_toml(joinpath(outdir, "fitted.toml"), cfg, scheme, res.x, total)
	write_convergence(joinpath(outdir, "convergence.csv"), res.history)
	for d in datasets
		write_curves(joinpath(outdir, "fit_$(safe_name(d.name)).csv"), scheme, d, kin)
	end
	isfile(checkpoint) && rm(checkpoint)
	open(joinpath(outdir, "summary.txt"), "w") do f
		write_summary(f, cfg, scheme, datasets, fp, res, per, total, opt)
	end
	write_summary(io, cfg, scheme, datasets, fp, res, per, total, opt)
	println(io, "Results written to ", outdir)
	return (x = res.x, names = scheme.names, cne = total, cne_per_dataset = per, result = res)
end

function write_summary(io, cfg, scheme, datasets, fp, res, per, total, opt)
	println(io, "TGAKinetics fit")
	println(io, "Run file: ", cfg["_path"])
	println(io, "Finished: ", Dates.format(now(), "yyyy-mm-dd HH:MM"))
	@printf(io, "Population %d, seed %d, %d generations, %d evaluations, %s on %d threads\n",
		opt.population, res.seed, res.generations, res.fevals, format_duration(res.seconds),
		Threads.nthreads())
	println(io, "Stopped: ", res.stop)
	opt.runs > 1 && println(io, "Best of $(opt.runs) runs: run $(res.run) (seed $(res.seed)); all runs in runs.csv")
	@printf(io, "\nCNE (tight solver tolerances): %.6e\n", total)
	for (d, w, c) in zip(datasets, fp.weights, per)
		@printf(io, "  %-20s weight %.4f  CNE %.6e\n", d.name, w, c)
	end
	println(io, "\nParameters (log10A in log10(1/s), E in kJ/mol):")
	for (i, nm) in enumerate(scheme.names)
		flag = near_bound(res.x[i], scheme.lower[i], scheme.upper[i]) ? "  <- at a bound" : ""
		@printf(io, "  %-28s %14.6g   [%g, %g]%s\n", nm, res.x[i], scheme.lower[i], scheme.upper[i], flag)
	end
end

"""
    simulate(path)

Run a scheme whose parameters are all fixed (for example `fitted.toml` from
`fit`) on every dataset: data files or set heating programs. Writes
`sim_<name>.csv` and, for datasets with data, reports the CNE.
"""
function simulate(path::AbstractString; io::IO = stdout)
	cfg, scheme, datasets = prepare(path)
	n_free(scheme) == 0 || cfgerror("simulate needs fixed values for every parameter; " *
									"these are ranges: $(join(scheme.names, ", "))")
	outdir = output_dir(cfg)
	N, M = length(scheme.components), length(scheme.reactions)
	kin = build_kinetics(scheme, Float64[], Val(N), Val(M))
	kin === nothing && cfgerror("the fixed yields or initial fractions are not feasible")
	obj = objective_settings(cfg)
	wm, wr = as_float(obj["mass_weight"], "mass_weight"), as_float(obj["mlr_weight"], "mlr_weight")
	wm, wr = wm / (wm + wr), wr / (wm + wr)
	describe(io, scheme)
	for d in datasets
		write_curves(joinpath(outdir, "sim_$(safe_name(d.name)).csv"), scheme, d, kin)
		if has_data(d)
			@printf(io, "%-20s CNE %.6e\n", d.name, cne_tight(kin, d, wm, wr))
		else
			println(io, d.name)
		end
	end
	println(io, "Curves written to ", outdir)
	return nothing
end
