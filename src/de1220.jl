# Self-adaptive differential evolution (the "DE 1220" algorithm of pagmo,
# Biscani and Izzo), as used in the original kinetics scripts: each member
# keeps its own F, CR and mutation variant, which are inherited with
# probability 0.9 and redrawn otherwise. One generation builds and scores all
# trial vectors in parallel from the previous population, then selects.
# Minimizes `f`.

const EXP_VARIANTS = (1, 2, 3, 4, 5, 11, 13, 15, 17)

struct DEOptions
	population::Int
	generations::Int
	runs::Int
	seed::Int
	variants::Vector{Int}
	ftol::Float64
	xtol::Float64
	report_every::Int
end

function DEOptions(d::AbstractDict)
	v = d["variants"]
	(v isa AbstractVector && !isempty(v) && all(x -> x isa Integer && 1 <= x <= 18, v)) ||
		cfgerror("optimizer.variants must be a list of integers from 1 to 18")
	o = DEOptions(as_int(d["population"], "population"), as_int(d["generations"], "generations"),
		as_int(d["runs"], "runs"), as_int(d["seed"], "seed"), Int.(v), as_float(d["ftol"], "ftol"),
		as_float(d["xtol"], "xtol"), as_int(d["report_every"], "report_every"))
	o.population >= 8 || cfgerror("optimizer.population must be at least 8")
	o.runs >= 1 || cfgerror("optimizer.runs must be at least 1")
	return o
end

# Element n of the mutant vector for `variant` (exponential and binomial
# versions of a strategy share the same formula).
@inline function mutant(variant, n, i, pop, best, F, r)
	x = pop
	v = variant
	if v == 1 || v == 6        # best/1
		return best[n] + F * (x[r[2]][n] - x[r[3]][n])
	elseif v == 2 || v == 7    # rand/1
		return x[r[1]][n] + F * (x[r[2]][n] - x[r[3]][n])
	elseif v == 3 || v == 8    # rand-to-best/1
		xi = x[i][n]
		return xi + F * (best[n] - xi) + F * (x[r[1]][n] - x[r[2]][n])
	elseif v == 4 || v == 9    # best/2
		return best[n] + F * (x[r[1]][n] - x[r[2]][n]) + F * (x[r[3]][n] - x[r[4]][n])
	elseif v == 5 || v == 10   # rand/2
		return x[r[5]][n] + F * (x[r[1]][n] - x[r[2]][n]) + F * (x[r[3]][n] - x[r[4]][n])
	elseif v == 11 || v == 12  # rand/3
		return x[r[1]][n] + F * (x[r[2]][n] - x[r[3]][n]) + F * (x[r[4]][n] - x[r[5]][n]) +
			   F * (x[r[6]][n] - x[r[7]][n])
	elseif v == 13 || v == 14  # best/3
		return best[n] + F * (x[r[2]][n] - x[r[3]][n]) + F * (x[r[4]][n] - x[r[5]][n]) +
			   F * (x[r[6]][n] - x[r[7]][n])
	elseif v == 15 || v == 16  # rand-to-current/2
		return x[r[1]][n] + F * (x[r[2]][n] - x[i][n]) + F * (x[r[3]][n] - x[r[4]][n])
	else                       # rand-to-best-and-current/2
		return x[r[1]][n] + F * (x[r[2]][n] - x[i][n]) - F * (x[r[3]][n] - best[n])
	end
end

function make_trial!(trial, rng, i, pop, best, Fi, CRi, Vi, opts, lower, upper, r, idxs)
	NP = length(pop)
	dim = length(trial)
	# seven distinct members other than i (Durstenfeld shuffle)
	for k in 1:NP
		idxs[k] = k
	end
	idxs[i], idxs[NP] = idxs[NP], idxs[i]
	for j in 1:7
		k = rand(rng, 1:(NP-j))
		r[j] = idxs[k]
		idxs[k], idxs[NP-j] = idxs[NP-j], idxs[k]
	end
	V = rand(rng) < 0.9 ? Vi : opts.variants[rand(rng, 1:length(opts.variants))]
	F = rand(rng) < 0.9 ? Fi : rand(rng) * 0.9 + 0.1
	CR = rand(rng) < 0.9 ? CRi : rand(rng)

	copyto!(trial, pop[i])
	n = rand(rng, 1:dim)
	if V in EXP_VARIANTS
		L = 0
		while true
			trial[n] = mutant(V, n, i, pop, best, F, r)
			n = mod(n, dim) + 1
			L += 1
			(rand(rng) < CR && L < dim) || break
		end
	else
		for L in 1:dim
			if rand(rng) < CR || L == dim
				trial[n] = mutant(V, n, i, pop, best, F, r)
			end
			n = mod(n, dim) + 1
		end
	end
	for j in 1:dim
		if trial[j] < lower[j] || trial[j] > upper[j]
			trial[j] = lower[j] + rand(rng) * (upper[j] - lower[j])
		end
	end
	return F, CR, V
end

"""
    sobol_population(lower, upper, NP, run=1)

Initial population from a Sobol sequence, skipping the first run × 2^m
points (2^m ≥ NP; Joe and Kuo 2003), so that repeat runs start from
different, equally well spread points.
"""
function sobol_population(lower, upper, NP, run = 1)
	dim = length(lower)
	seq = SobolSeq(dim)
	Sobol.skip(seq, run * 2^ceil(Int, log2(max(NP, 2))); exact = true)
	return [lower .+ next!(seq) .* (upper .- lower) for _ in 1:NP]
end

"""
    de1220(f, lower, upper, opts; run=1, on_report=nothing) -> NamedTuple

Minimize `f` within the box `[lower, upper]`. Repeat run `run` uses seed
`opts.seed + run - 1` and its own Sobol start. `on_report(gen, best_x,
best_f)` is called at every report line (used to save a checkpoint).
"""
function de1220(f, lower::Vector{Float64}, upper::Vector{Float64}, opts::DEOptions;
	run::Int = 1, on_report = nothing, io::IO = stdout)
	NP = opts.population
	dim = length(lower)
	seed = opts.seed + run - 1
	rng = Xoshiro(seed)
	t_start = time()

	pop = sobol_population(lower, upper, NP, run)
	fit = Vector{Float64}(undef, NP)
	Threads.@threads :greedy for i in 1:NP
		fit[i] = f(pop[i])
	end
	fevals = NP
	Fm = [rand(rng) * 0.9 + 0.1 for _ in 1:NP]
	CRm = [rand(rng) for _ in 1:NP]
	Vm = [opts.variants[rand(rng, 1:length(opts.variants))] for _ in 1:NP]

	ib = argmin(fit)
	best_x, best_f = copy(pop[ib]), fit[ib]
	history = [(0, fevals, best_f, time() - t_start)]
	@printf(io, "Initial population: %d members, best CNE %.6e\n", NP, best_f)

	trials = [Vector{Float64}(undef, dim) for _ in 1:NP]
	tF = zeros(NP)
	tCR = zeros(NP)
	tV = zeros(Int, NP)
	tfit = zeros(NP)
	rs = [Vector{Int}(undef, 7) for _ in 1:NP]
	idxs = [Vector{Int}(undef, NP) for _ in 1:NP]
	rngs = [Xoshiro(0) for _ in 1:NP]
	stop = "generation limit reached"
	gen = 0
	gen_times = Float64[]

	for g in 1:opts.generations
		gen = g
		tg = time()
		gbest = copy(best_x)
		Threads.@threads :greedy for i in 1:NP
			Random.seed!(rngs[i], hash((seed, g, i)))
			tF[i], tCR[i], tV[i] = make_trial!(trials[i], rngs[i], i, pop, gbest, Fm[i], CRm[i], Vm[i],
				opts, lower, upper, rs[i], idxs[i])
			tfit[i] = f(trials[i])
		end
		fevals += NP
		for i in 1:NP
			if tfit[i] <= fit[i] || isnan(fit[i])
				copyto!(pop[i], trials[i])
				fit[i] = tfit[i]
				Fm[i], CRm[i], Vm[i] = tF[i], tCR[i], tV[i]
				if tfit[i] <= best_f
					best_f = tfit[i]
					copyto!(best_x, trials[i])
				end
			end
		end
		push!(gen_times, time() - tg)
		push!(history, (g, fevals, best_f, time() - t_start))

		ib, iw = argmin(fit), argmax(fit)
		dx = sum(abs, pop[iw] .- pop[ib])
		df = abs(fit[iw] - fit[ib])
		if opts.report_every > 0 && (g % opts.report_every == 0 || g == 1)
			recent = mean(@view gen_times[max(1, end - opts.report_every + 1):end])
			@printf(io, "gen %6d  evals %8d  best CNE %.6e  spread %.2e  %.2f s/gen  left %s\n",
				g, fevals, best_f, df, recent, format_duration(recent * (opts.generations - g)))
			flush(io)
			on_report === nothing || on_report(g, best_x, best_f)
		end
		if dx < opts.xtol
			stop = "population collapsed in parameter space (xtol)"
			break
		elseif df < opts.ftol
			stop = "population collapsed in CNE (ftol)"
			break
		end
	end
	return (x = best_x, f = best_f, generations = gen, fevals = fevals, stop = stop,
		history = history, seconds = time() - t_start, run = run, seed = seed)
end

function format_duration(s)
	s = round(Int, s)
	h, r = divrem(s, 3600)
	m, sec = divrem(r, 60)
	return @sprintf("%d:%02d:%02d", h, m, sec)
end
