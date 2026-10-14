# Reading the TOML run file. Unknown keys are errors so that typos do not
# silently fall back to defaults.

const TOP_KEYS = Set(["output_dir", "components", "initial", "reaction", "data",
	"dataset", "objective", "optimizer", "solver"])
const DATA_KEYS = Set(["header_row", "time_column", "temperature_column", "mass_column",
	"mlr_column", "time_unit", "temperature_unit", "reference", "dry_window",
	"fit_range", "grid_step", "smoothing_segments", "smoothing_lambda"])
const DATASET_KEYS = union(DATA_KEYS, Set(["name", "file", "weight",
	"heating_rate", "T_start", "T_end"]))
const REACTION_KEYS = Set(["name", "reactant", "products", "log10A", "E", "order"])
const OBJECTIVE_KEYS = Set(["mlr_weight", "mass_weight"])
const OPTIMIZER_KEYS = Set(["population", "generations", "runs", "seed", "report_every",
	"ftol", "xtol", "variants"])
const SOLVER_KEYS = Set(["abstol", "reltol", "maxiters"])

const DATA_DEFAULTS = Dict{String, Any}(
	"header_row" => 1,
	"time_column" => "time",
	"temperature_column" => "temperature",
	"mass_column" => "mass",
	"mlr_column" => "",
	"time_unit" => "s",
	"temperature_unit" => "K",
	"reference" => "dry",
	"dry_window" => [373.0, 423.0],
	"fit_range" => [0.0, 1.0e4],
	"grid_step" => 1.0,
	"smoothing_segments" => 500,
	"smoothing_lambda" => "auto",
)

const OBJECTIVE_DEFAULTS = Dict{String, Any}("mlr_weight" => 0.6, "mass_weight" => 0.4)

const OPTIMIZER_DEFAULTS = Dict{String, Any}(
	"population" => 50,
	"generations" => 1000,
	"runs" => 1,
	"seed" => 1,
	"report_every" => 10,
	"ftol" => 1e-8,
	"xtol" => 1e-8,
	"variants" => [2, 3, 7, 10, 13, 14, 15, 16],
)

const SOLVER_DEFAULTS = Dict{String, Any}("abstol" => 1e-6, "reltol" => 1e-4,
	"maxiters" => 100_000)

struct ConfigError <: Exception
	msg::String
end
Base.showerror(io::IO, e::ConfigError) = print(io, "Config error: ", e.msg)
cfgerror(msg) = throw(ConfigError(msg))

function check_keys(d::AbstractDict, allowed::Set{String}, where::AbstractString)
	for k in keys(d)
		k in allowed || cfgerror("unknown key \"$k\" in $where. Allowed: " *
								 join(sort(collect(allowed)), ", "))
	end
end

"""
    load_config(path) -> Dict

Read a run file and return it as a Dict, with `"_dir"` set to the folder that
holds the file (relative paths in the file are taken from there).
"""
function load_config(path::AbstractString)
	isfile(path) || cfgerror("file not found: $(abspath(path))")
	cfg = try
		TOML.parsefile(path)
	catch e
		cfgerror("could not parse $(abspath(path)) as TOML: $(sprint(showerror, e))")
	end
	check_keys(cfg, TOP_KEYS, "the top level of the run file")
	cfg["_dir"] = dirname(abspath(path))
	cfg["_path"] = abspath(path)
	return cfg
end

resolve_path(cfg, p::AbstractString) = isabspath(p) ? p : normpath(joinpath(cfg["_dir"], p))

function section(cfg, name, allowed, defaults)
	d = get(cfg, name, Dict{String, Any}())
	d isa AbstractDict || cfgerror("[$name] must be a table")
	check_keys(d, allowed, "[$name]")
	return merge(defaults, d)
end

objective_settings(cfg) = section(cfg, "objective", OBJECTIVE_KEYS, OBJECTIVE_DEFAULTS)
optimizer_settings(cfg) = section(cfg, "optimizer", OPTIMIZER_KEYS, OPTIMIZER_DEFAULTS)
solver_settings(cfg) = section(cfg, "solver", SOLVER_KEYS, SOLVER_DEFAULTS)

function output_dir(cfg)
	dir = resolve_path(cfg, string(get(cfg, "output_dir", "results")))
	mkpath(dir)
	return dir
end

as_float(x, what) = x isa Real ? Float64(x) : cfgerror("$what must be a number, got $(repr(x))")
as_int(x, what) = x isa Integer ? Int(x) : cfgerror("$what must be an integer, got $(repr(x))")

function as_range(x, what)
	(x isa AbstractVector && length(x) == 2 && all(v -> v isa Real, x)) ||
		cfgerror("$what must be [lower, upper], got $(repr(x))")
	lo, hi = Float64(x[1]), Float64(x[2])
	lo < hi || cfgerror("$what: lower ($lo) must be below upper ($hi)")
	return lo, hi
end
