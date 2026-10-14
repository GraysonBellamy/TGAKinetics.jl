using TGAKinetics
using Test
using TOML
using StaticArrays
using Random
using Printf

const TK = TGAKinetics
const EXAMPLE = joinpath(@__DIR__, "..", "examples", "synthetic")

scheme_from(str) = TK.parse_scheme(TOML.parse(str))

function kin_of(s::TK.Scheme, x = Float64[])
	N, M = length(s.components), length(s.reactions)
	return TK.build_kinetics(s, x, Val(N), Val(M))
end

constant_T(T, tend) = TK.TempProgram(collect(range(0.0, tend; length = 101)), fill(T, 101); nseg = 10, lambda = 1.0)
ramp_T(T0, beta, tend) = (tl = collect(range(0.0, tend; length = 101));
TK.TempProgram(tl, T0 .+ beta .* tl; nseg = 10, lambda = 1.0))

@testset "P-spline smoothing" begin
	rng = Xoshiro(1)
	t = collect(0.0:0.5:600.0)
	f(x) = 1 - 0.5 / (1 + exp(-(x - 300) / 30))
	df(x) = -0.5 * exp(-(x - 300) / 30) / 30 / (1 + exp(-(x - 300) / 30))^2
	y = f.(t) .+ 1e-4 .* randn(rng, length(t))
	ps = TK.pspline_fit(t, y)
	inner = 50:(length(t)-50)
	@test maximum(abs.(TK.pspline_eval(ps, t)[inner] .- f.(t)[inner])) < 2e-4
	@test maximum(abs.(TK.pspline_eval(ps, t; derivative = true)[inner] .- df.(t)[inner])) < 0.03 * maximum(abs, df.(t))
	# a straight line is reproduced exactly for any smoothing weight
	line = TK.pspline_fit(t, 300 .+ 0.2 .* t; nseg = 10, lambda = 1e3)
	@test maximum(abs.(TK.pspline_eval(line, t) .- (300 .+ 0.2 .* t))) < 1e-8
	@test all(abs.(TK.pspline_eval(line, t; derivative = true) .- 0.2) .< 1e-9)
end

@testset "mass-loss rate differences" begin
	t = sort!(rand(Xoshiro(2), 40)) .* 100
	m = 3 .+ 0.2 .* t .- 0.01 .* t .^ 2       # second-order differences are exact on a quadratic
	mlr = TK.mass_loss_rate!(similar(m), m, t)
	@test maximum(abs.(mlr .- (-(0.2 .- 0.02 .* t)))) < 1e-9
end

@testset "exact solutions" begin
	s = scheme_from("""
	components = ["a", "b"]
	[[reaction]]
	reactant = "a"
	products = { b = 0.3 }
	log10A = 5.0
	E = 80.0
	""")
	kin = kin_of(s)
	k = 1e5 * exp(-80e3 / (TK.R_GAS * 600.0))
	ts = collect(range(0.0, 3 / k; length = 50))
	u, ok = TK.simulate_run(kin, constant_T(600.0, ts[end]), ts, TK.TIGHT)
	@test ok
	@test maximum(abs(v[1] - exp(-k * t)) for (v, t) in zip(u, ts)) < 1e-7
	@test maximum(abs(sum(v) - (0.3 + 0.7 * exp(-k * t))) for (v, t) in zip(u, ts)) < 1e-7

	# constant heating rate: ln a = -∫ A exp(-E / R T(t)) dt
	beta = 10 / 60
	tend = 400 / beta
	tf = range(0.0, tend; length = 400_001)
	ts = collect(tf[1:2000:end])
	u, ok = TK.simulate_run(kin, ramp_T(450.0, beta, tend), ts, TK.TIGHT)
	@test ok
	rate = 1e5 .* exp.(-80e3 ./ (TK.R_GAS .* (450.0 .+ beta .* tf)))
	cum = [0.0; cumsum((rate[1:end-1] .+ rate[2:end]) ./ 2 .* step(tf))]
	exact = exp.(-cum[1:2000:end])
	@test maximum(abs.(first.(u) .- exact)) < 1e-6

	# second order: da/dt = -k a^2  =>  a = 1 / (1 + k t)
	s2 = scheme_from("""
	components = ["a"]
	[[reaction]]
	reactant = "a"
	log10A = 5.0
	E = 80.0
	order = 2.0
	""")
	ts = collect(range(0.0, 10 / k; length = 50))
	u, ok = TK.simulate_run(kin_of(s2), constant_T(600.0, ts[end]), ts, TK.TIGHT)
	@test ok
	@test maximum(abs(v[1] - 1 / (1 + k * t)) for (v, t) in zip(u, ts)) < 1e-7

	# parallel components with initial fractions and a "rest" entry
	s3 = scheme_from("""
	components = ["p", "q", "c"]
	[initial]
	p = 0.4
	q = "rest"
	[[reaction]]
	reactant = "p"
	log10A = 5.0
	E = 80.0
	[[reaction]]
	reactant = "q"
	products = { c = 0.25 }
	log10A = 6.0
	E = 80.0
	""")
	ts = collect(range(0.0, 3 / k; length = 50))
	u, ok = TK.simulate_run(kin_of(s3), constant_T(600.0, ts[end]), ts, TK.TIGHT)
	@test ok
	exact3 = [0.4 * exp(-k * t) + 0.6 * (0.25 + 0.75 * exp(-10k * t)) for t in ts]
	@test maximum(abs.(sum.(u) .- exact3)) < 1e-7
end

@testset "Jacobian and time derivative" begin
	s = scheme_from("""
	components = ["a", "b", "c", "d"]
	[[reaction]]
	reactant = "a"
	products = { b = 0.6, c = 0.1 }
	log10A = 8.0
	E = 110.0
	order = 1.5
	[[reaction]]
	reactant = "b"
	products = { d = 0.4 }
	log10A = 10.0
	E = 140.0
	[[reaction]]
	reactant = "c"
	log10A = 3.0
	E = 60.0
	order = 0.7
	""")
	kin = kin_of(s)
	p = TK.ODEParams(kin, ramp_T(500.0, 0.2, 2000.0))
	u = SVector(0.5, 0.2, 0.1, 0.05)
	t = 700.0
	J = TK.jac(u, p, t)
	h = 1e-7
	Jfd = hcat([(TK.rhs(u + h * SVector(ntuple(i -> i == a ? 1.0 : 0.0, 4)), p, t) .-
				 TK.rhs(u - h * SVector(ntuple(i -> i == a ? 1.0 : 0.0, 4)), p, t)) ./ 2h for a in 1:4]...)
	@test maximum(abs.(J .- Jfd)) < 1e-6 * max(1, maximum(abs, J))
	tg = TK.tgrad(u, p, t)
	tgfd = (TK.rhs(u, p, t + 1e-3) .- TK.rhs(u, p, t - 1e-3)) ./ 2e-3
	@test maximum(abs.(tg .- tgfd)) < 1e-6 * max(1e-12, maximum(abs, tg))
end

@testset "config checks" begin
	@test_throws TK.ConfigError scheme_from("""
	components = ["a", "b"]
	[[reaction]]
	reactant = "a"
	products = { b = 0.6 }
	log10A = 5.0
	E = 80.0
	Ea = 3.0
	""")
	@test_throws TK.ConfigError scheme_from("""
	components = ["a", "b", "c"]
	[[reaction]]
	reactant = "a"
	products = { b = 0.6, c = 0.5 }
	log10A = 5.0
	E = 80.0
	""")
	@test_throws TK.ConfigError scheme_from("""
	components = ["a", "b"]
	[initial]
	a = [0.2, 0.5]
	b = 0.5
	[[reaction]]
	reactant = "a"
	log10A = 5.0
	E = 80.0
	""")
	@test_throws TK.ConfigError scheme_from("""
	components = ["a"]
	[[reaction]]
	reactant = "x"
	log10A = 5.0
	E = 80.0
	""")
	s = scheme_from("""
	components = ["a", "b", "c"]
	[[reaction]]
	reactant = "a"
	products = { b = [0.5, 0.9], c = [0.2, 0.4] }
	log10A = [1.0, 10.0]
	E = 80.0
	""")
	@test s.names == ["R1.log10A", "R1.yield.b", "R1.yield.c"]
	# yields that add up to more than 1 are rejected at evaluation time
	@test kin_of(s, [5.0, 0.9, 0.4]) === nothing
	@test kin_of(s, [5.0, 0.6, 0.4]) !== nothing
end

@testset "example: check, fit, simulate" begin
	dir = mktempdir()
	cfgtext = read(joinpath(EXAMPLE, "run.toml"), String)
	cfgtext = replace(cfgtext, "output_dir = \"results\"" => "output_dir = \"$(escape_string(dir))\"",
		"file = \"data/" => "file = \"$(escape_string(joinpath(EXAMPLE, "data")))/",
		"population = 50" => "population = 16", "generations = 2000" => "generations = 20", "runs = 4" => "runs = 2",
		"report_every = 500" => "report_every = 0")
	run = joinpath(dir, "run.toml")
	write(run, cfgtext)

	cfg, scheme, ds = TK.prepare(run)
	@test [d.name for d in ds] == ["5Kmin", "10Kmin", "20Kmin"]
	for d in ds
		@test isapprox(d.ref_mass, 5.0 * 0.94; rtol = 2e-3)   # dry mass of the synthetic sample
	end
	fp = TK.FitProblem(scheme, ds, TK.objective_settings(cfg), TK.SolverOptions(TK.solver_settings(cfg)))
	truth = [9.17318626841, 119.2, 0.81, 12.296, 160.6, 0.82, 13.6985, 189.6, 0.39,
		4.62, 92.26, 0.81, 0.587, 52.57, 0.80]
	@test TK.objective(fp, truth) < 1e-3

	TK.check(run; io = devnull)
	@test isfile(joinpath(dir, "data_10Kmin.csv"))

	r1 = TK.fit(run; io = devnull)
	r2 = TK.fit(run; io = devnull)
	@test r1.x == r2.x                      # same seed, same answer
	for f in ("parameters.csv", "fitted.toml", "convergence.csv", "summary.txt", "fit_5Kmin.csv")
		@test isfile(joinpath(dir, f))
	end
	TK.simulate(joinpath(dir, "fitted.toml"); io = devnull)
	sim = TK.CSV.File(joinpath(dir, "simulation", "sim_10Kmin.csv"))
	fitc = TK.CSV.File(joinpath(dir, "fit_10Kmin.csv"))
	@test maximum(abs.(sim.mass_model .- fitc.mass_model)) < 1e-9
end
