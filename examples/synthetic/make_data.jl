# Make the synthetic TGA files in data/ from true_scheme.toml: 5 mg samples,
# one reading per second, balance noise 1 µg, thermocouple noise 0.05 K.
#
#   julia --project=../.. make_data.jl

using TGAKinetics, Random, Printf
const T = TGAKinetics

cfg = T.load_config(joinpath(@__DIR__, "true_scheme.toml"))
scheme = T.parse_scheme(cfg)
N, M = length(scheme.components), length(scheme.reactions)
kin = T.build_kinetics(scheme, Float64[], Val(N), Val(M))
rng = Xoshiro(2026)
m0 = 5.0  # mg

for (k, entry) in enumerate(cfg["dataset"])
	d = T.load_dataset(cfg, entry, k)
	ts = collect(0.0:1.0:d.t[end])
	u, ok = T.simulate_run(kin, d.prog, ts, T.TIGHT)
	ok || error("solve failed")
	open(joinpath(@__DIR__, "data", "tga_$(d.name).csv"), "w") do io
		println(io, "Time (min),Temperature (C),Mass (mg)")
		for (t, v) in zip(ts, u)
			Tc = T.temperature(d.prog, t) - 273.15 + 0.05 * randn(rng)
			m = m0 * sum(v) + 0.001 * randn(rng)
			@printf(io, "%.5f,%.3f,%.5f\n", t / 60, Tc, m)
		end
	end
	println("wrote tga_$(d.name).csv ($(length(ts)) rows)")
end
