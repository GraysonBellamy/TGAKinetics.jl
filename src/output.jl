# Result files.

safe_name(s) = replace(s, r"[^A-Za-z0-9_.-]+" => "_")

function write_table(path, names::Vector{String}, cols::Vector)
	nt = NamedTuple{Tuple(Symbol.(names))}(Tuple(cols))
	CSV.write(path, nt)
end

"""
Measured and model curves at every compared point of dataset `d`.
"""
function write_curves(path, scheme::Scheme, d::Dataset, kin::Kinetics)
	mass, mlr, comps = model_curves(kin, d)
	names = ["time_s", "temperature_K", "mass_exp", "mass_model", "mlr_exp_per_s",
		"mlr_model_per_s", "in_fit_range"]
	cols = Any[d.t, d.T, d.mass, mass, d.mlr, mlr, Int.(d.infit)]
	for (c, v) in zip(scheme.components, comps)
		push!(names, "comp_" * c)
		push!(cols, v)
	end
	write_table(path, names, cols)
end

function write_processed(path, d::Dataset)
	p = d.processed
	write_table(path, ["time_s", "temperature_K", "temperature_smoothed_K", "mass", "mass_smoothed", "mlr_per_s"],
		Any[p.time_s, p.temperature_K, p.temperature_smoothed_K, p.mass, p.mass_smoothed, p.mlr_per_s])
end

near_bound(v, lo, hi) = min(v - lo, hi - v) < 0.005 * (hi - lo)

function write_parameters(path, scheme::Scheme, x)
	open(path, "w") do io
		println(io, "parameter,value,lower,upper,near_bound")
		for (i, n) in enumerate(scheme.names)
			lo, hi = scheme.lower[i], scheme.upper[i]
			@printf(io, "%s,%.10g,%.10g,%.10g,%d\n", n, x[i], lo, hi, near_bound(x[i], lo, hi))
		end
	end
end

function write_convergence(path, history)
	open(path, "w") do io
		println(io, "generation,evaluations,best_cne,elapsed_s")
		for (g, n, f, t) in history
			@printf(io, "%d,%d,%.10e,%.3f\n", g, n, f, t)
		end
	end
end

function write_runs(path, results)
	open(path, "w") do io
		println(io, "run,seed,best_cne,generations,evaluations,seconds,stop")
		for r in results
			@printf(io, "%d,%d,%.10e,%d,%d,%.1f,\"%s\"\n", r.run, r.seed, r.f, r.generations,
				r.fevals, r.seconds, r.stop)
		end
	end
end

fmt(v::Float64) = repr(v)

"""
Write a complete run file with the fitted values filled in, so that
`simulate` can be run on it directly.
"""
function write_fitted_toml(path, cfg, scheme::Scheme, x, total_cne)
	open(path, "w") do io
		println(io, "# Written by TGAKinetics on $(Dates.format(now(), "yyyy-mm-dd HH:MM")) from")
		println(io, "# $(cfg["_path"])")
		@printf(io, "# CNE = %.6e (tight solver tolerances)\n\n", total_cne)
		println(io, "output_dir = \"simulation\"")
		println(io, "components = [", join(("\"$c\"" for c in scheme.components), ", "), "]\n")
		println(io, "[initial]")
		for (i, c) in enumerate(scheme.components)
			if i == scheme.rest
				println(io, "$c = \"rest\"")
			else
				v = slotvalue(scheme.initial[i], x)
				v != 0 && println(io, "$c = ", fmt(v))
			end
		end
		for r in scheme.reactions
			println(io, "\n[[reaction]]")
			println(io, "name = \"$(r.name)\"")
			println(io, "reactant = \"$(scheme.components[r.reactant])\"")
			prods = join(("$(scheme.components[p]) = $(fmt(slotvalue(y, x)))" for (p, y) in zip(r.products, r.yields)), ", ")
			println(io, isempty(prods) ? "products = {}" : "products = { $prods }")
			println(io, "log10A = ", fmt(slotvalue(r.log10A, x)))
			println(io, "E = ", fmt(slotvalue(r.E, x)))
			println(io, "order = ", fmt(slotvalue(r.order, x)))
		end
		println(io)
		haskey(cfg, "data") && TOML.print(io, Dict("data" => cfg["data"]))
		sets = map(cfg["dataset"]) do e
			e = copy(e)
			haskey(e, "file") && (e["file"] = resolve_path(cfg, e["file"]))
			e
		end
		println(io)
		TOML.print(io, Dict("dataset" => sets))
		for key in ("objective", "solver")
			if haskey(cfg, key)
				println(io)
				TOML.print(io, Dict(key => cfg[key]))
			end
		end
	end
end
