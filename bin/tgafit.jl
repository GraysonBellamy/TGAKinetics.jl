# Command-line entry point.
#
#   julia -t auto bin/tgafit.jl check    run.toml
#   julia -t auto bin/tgafit.jl fit      run.toml
#   julia -t auto bin/tgafit.jl simulate run.toml

Base.set_active_project(joinpath(@__DIR__, "..", "Project.toml"))

using TGAKinetics

const USAGE = """
usage: julia -t auto bin/tgafit.jl <command> <run file>

commands:
  check     read the data and scheme, write the processed data, estimate run time
  fit       fit the free parameters
  simulate  run a scheme with fixed parameters (e.g. results/fitted.toml)
"""

function main(args)
	if length(args) != 2 || !(args[1] in ("check", "fit", "simulate"))
		print(stderr, USAGE)
		return 1
	end
	cmd, path = args
	try
		cmd == "check" ? check(path) : cmd == "fit" ? fit(path) : simulate(path)
	catch e
		e isa TGAKinetics.ConfigError || rethrow()
		println(stderr, sprint(showerror, e))
		return 2
	end
	return 0
end

exit(main(ARGS))
