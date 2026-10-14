# Loading one TGA run: read the columns, convert units, get the mass-loss rate,
# normalize by the reference mass, and set the points the objective compares.

"""
Sample temperature against time as a smooth cubic spline (the P-spline fit
of the measured temperature), so the model sees no kinks from thermocouple
noise. Outside the measured span the temperature is held at the end values.
"""
struct TempProgram
	t0::Float64
	t1::Float64
	spline::PSplineFit
end

@inline function temperature(p::TempProgram, t)
	tc = clamp(t, p.t0, p.t1)
	return spline_value(p.spline, tc)
end

@inline function temperature_slope(p::TempProgram, t)
	(t < p.t0 || t > p.t1) && return 0.0
	return spline_slope(p.spline, t)
end

function TempProgram(t::AbstractVector, T::AbstractVector; nseg::Int = 500, lambda = :auto)
	return TempProgram(t[1], t[end], pspline_fit(t, T; nseg = nseg, lambda = lambda))
end

struct Dataset
	name::String
	file::String              # "" for a set heating program
	weight::Float64
	prog::TempProgram
	t::Vector{Float64}        # compared points, s (all of them, for output)
	T::Vector{Float64}        # K
	mass::Vector{Float64}     # measured, normalized (NaN without data)
	mlr::Vector{Float64}      # measured, 1/s
	infit::BitVector          # inside fit_range
	mass_max::Float64
	mlr_max::Float64
	ref_mass::Float64         # in the file's mass unit
	ref_temperature::Float64  # K
	lambda::Float64           # P-spline smoothing weight (NaN if not used)
	processed::Any            # full-resolution series for data_<name>.csv
end

has_data(d::Dataset) = !isempty(d.file)

function linear_interp(x::AbstractVector, y::AbstractVector, xq::AbstractVector)
	out = Vector{Float64}(undef, length(xq))
	j = 1
	n = length(x)
	@inbounds for (k, q) in enumerate(xq)
		if q <= x[1]
			out[k] = y[1]
		elseif q >= x[n]
			out[k] = y[n]
		else
			while x[j+1] < q
				j += 1
			end
			while j > 1 && x[j] > q
				j -= 1
			end
			f = (q - x[j]) / (x[j+1] - x[j])
			out[k] = (1 - f) * y[j] + f * y[j+1]
		end
	end
	return out
end

function dataset_settings(cfg, entry, k)
	check_keys(entry, DATASET_KEYS, "[[dataset]] number $k")
	defaults = get(cfg, "data", Dict{String, Any}())
	defaults isa AbstractDict || cfgerror("[data] must be a table")
	check_keys(defaults, DATA_KEYS, "[data]")
	return merge(DATA_DEFAULTS, defaults, entry)
end

function read_column(file, csv, name, what)
	name in string.(propertynames(csv)) ||
		cfgerror("$(basename(file)) has no $what column \"$name\". Columns: " *
				 join(string.(propertynames(csv)), ", "))
	col = getproperty(csv, Symbol(name))
	eltype(col) <: Union{Missing, Real} ||
		cfgerror("column \"$name\" of $(basename(file)) is not numeric. If there is a units row " *
				 "or text below the header, remove it; if the header is not on line 1, set header_row")
	return [ismissing(v) ? NaN : Float64(v) for v in col]
end

"""
    load_dataset(cfg, entry, k) -> Dataset

Process one `[[dataset]]` entry (a data file, or a set heating program for
`simulate`).
"""
function load_dataset(cfg, entry, k)
	s = dataset_settings(cfg, entry, k)
	name = string(get(s, "name", "run$k"))
	weight = as_float(get(s, "weight", 1.0), "weight of dataset $name")
	weight >= 0 || cfgerror("weight of dataset $name must be >= 0")
	grid_step = as_float(s["grid_step"], "grid_step")
	fit_lo, fit_hi = as_range(s["fit_range"], "fit_range")

	if !haskey(s, "file")
		return heating_program(name, s, weight, grid_step, fit_lo, fit_hi)
	end

	file = resolve_path(cfg, string(s["file"]))
	isfile(file) || cfgerror("data file not found: $file")
	hr = as_int(s["header_row"], "header_row")
	csv = CSV.File(file; normalizenames = false, header = hr, skipto = hr + 1)
	t = read_column(file, csv, s["time_column"], "time")
	T = read_column(file, csv, s["temperature_column"], "temperature")
	m = read_column(file, csv, s["mass_column"], "mass")
	has_mlr = !isempty(s["mlr_column"])
	r = has_mlr ? read_column(file, csv, s["mlr_column"], "mass-loss rate") : zeros(length(t))

	tu = s["time_unit"]
	tu in ("s", "min") || cfgerror("time_unit must be \"s\" or \"min\"")
	tu == "min" && (t .*= 60.0; r ./= 60.0)
	Tu = s["temperature_unit"]
	Tu in ("K", "C") || cfgerror("temperature_unit must be \"K\" or \"C\"")
	Tu == "C" && (T .+= 273.15)

	keep = isfinite.(t) .& isfinite.(T) .& isfinite.(m) .& isfinite.(r)
	t, T, m, r = t[keep], T[keep], m[keep], r[keep]
	order = sortperm(t)
	t, T, m, r = t[order], T[order], m[order], r[order]
	uniq = [true; diff(t) .> 0]
	t, T, m, r = t[uniq], T[uniq], m[uniq], r[uniq]
	length(t) >= 20 || cfgerror("$(basename(file)) has fewer than 20 usable rows")

	lam = s["smoothing_lambda"]
	lam = lam == "auto" ? :auto : as_float(lam, "smoothing_lambda")
	ps = pspline_fit(t, m; nseg = as_int(s["smoothing_segments"], "smoothing_segments"), lambda = lam)
	ms = pspline_eval(ps, t)
	has_mlr || (r = -pspline_eval(ps, t; derivative = true))

	ref = s["reference"]
	if ref == "initial"
		i0 = 1
	elseif ref == "dry"
		lo, hi = as_range(s["dry_window"], "dry_window")
		inwin = findall(x -> lo < x < hi, T)
		isempty(inwin) && cfgerror("$(basename(file)): no points between $lo and $hi K for the dry-mass reference")
		i0 = inwin[argmin(abs.(r[inwin]))]
	else
		cfgerror("reference must be \"dry\" or \"initial\", got $(repr(ref))")
	end
	ref_mass = ms[i0]
	ref_mass > 0 || cfgerror("$(basename(file)): reference mass is $ref_mass")

	t, T, m, ms, r = t[i0:end], T[i0:end], m[i0:end], ms[i0:end], r[i0:end]
	mn, msn, rn = m ./ ref_mass, ms ./ ref_mass, r ./ ref_mass

	prog = TempProgram(t, T; nseg = as_int(s["smoothing_segments"], "smoothing_segments"))
	span = t[end] - t[1]
	mean_rate = sum(abs, diff(T)) / span
	dt_eval = min(grid_step / max(mean_rate, eps()), span / 500)
	te = collect(range(t[1], t[end]; step = dt_eval))
	Te = [temperature(prog, x) for x in te]
	me = linear_interp(t, mn, te)
	re = linear_interp(t, rn, te)
	infit = BitVector(fit_lo .<= Te .<= fit_hi)
	count(infit) >= 10 || cfgerror("$name: fewer than 10 compared points inside fit_range [$fit_lo, $fit_hi] K")

	processed = (time_s = t, temperature_K = T, temperature_smoothed_K = [temperature(prog, x) for x in t],
		mass = mn, mass_smoothed = msn, mlr_per_s = rn)
	return Dataset(name, file, weight, prog, te, Te, me, re, infit,
		maximum(me[infit]), maximum(abs, re[infit]), ref_mass, T[1],
		has_mlr ? NaN : ps.lambda, processed)
end

function heating_program(name, s, weight, grid_step, fit_lo, fit_hi)
	for key in ("heating_rate", "T_start", "T_end")
		haskey(s, key) || cfgerror("dataset $name needs either `file` or heating_rate, T_start and T_end")
	end
	beta = as_float(s["heating_rate"], "heating_rate") / 60.0
	T0 = as_float(s["T_start"], "T_start")
	T1 = as_float(s["T_end"], "T_end")
	(beta > 0 && T1 > T0) || cfgerror("dataset $name: need heating_rate > 0 and T_end > T_start")
	tend = (T1 - T0) / beta
	tl = collect(range(0.0, tend; length = 101))
	prog = TempProgram(tl, T0 .+ beta .* tl; nseg = 10, lambda = 1.0)
	te = collect(range(0.0, tend; step = grid_step / beta))
	Te = T0 .+ beta .* te
	infit = BitVector(fit_lo .<= Te .<= fit_hi)
	nan = fill(NaN, length(te))
	return Dataset(name, "", weight, prog, te, Te, nan, copy(nan), infit, NaN, NaN, NaN, T0, NaN, nothing)
end
