# The reaction scheme: components, initial composition and reactions. Every
# scalar is either fixed (a number in the run file) or free ([lower, upper]).
# A `Slot` records which: `idx == 0` means fixed at `value`, otherwise the
# value is element `idx` of the optimizer's parameter vector.

struct Slot
	idx::Int
	value::Float64
end

@inline slotvalue(s::Slot, x) = s.idx == 0 ? s.value : @inbounds x[s.idx]

struct ReactionSpec
	name::String
	reactant::Int
	products::Vector{Int}
	log10A::Slot
	E::Slot        # kJ/mol
	order::Slot
	yields::Vector{Slot}
end

struct Scheme
	components::Vector{String}
	reactions::Vector{ReactionSpec}
	initial::Vector{Slot}   # one per component; ignored for `rest`
	rest::Int               # component that takes 1 - sum(others), or 0
	names::Vector{String}   # free parameters, in optimizer order
	lower::Vector{Float64}
	upper::Vector{Float64}
end

n_free(s::Scheme) = length(s.names)

mutable struct FreeList
	names::Vector{String}
	lower::Vector{Float64}
	upper::Vector{Float64}
end

function parse_slot!(fl::FreeList, x, name::String; lo = -Inf, hi = Inf)
	if x isa Real
		v = Float64(x)
		lo <= v <= hi || cfgerror("$name = $v is outside the allowed range [$lo, $hi]")
		return Slot(0, v)
	elseif x isa AbstractVector
		a, b = as_range(x, name)
		(lo <= a && b <= hi) || cfgerror("$name bounds [$a, $b] must lie within [$lo, $hi]")
		push!(fl.names, name)
		push!(fl.lower, a)
		push!(fl.upper, b)
		return Slot(length(fl.names), NaN)
	end
	cfgerror("$name must be a number (fixed) or [lower, upper] (fitted), got $(repr(x))")
end

function component_index(components, name, where)
	i = findfirst(==(name), components)
	i === nothing && cfgerror("$where refers to \"$name\", which is not in components = $(components)")
	return i
end

"""
    parse_scheme(cfg) -> Scheme

Build the reaction scheme from the `components`, `[initial]` and
`[[reaction]]` entries of a run file.
"""
function parse_scheme(cfg)
	haskey(cfg, "components") || cfgerror("missing `components = [...]` (put it above the first [table])")
	comps = cfg["components"]
	(comps isa AbstractVector && !isempty(comps) && all(c -> c isa AbstractString, comps)) ||
		cfgerror("components must be a list of names, e.g. [\"wood\", \"char\"]")
	comps = String.(comps)
	allunique(comps) || cfgerror("component names must be unique: $(comps)")

	fl = FreeList(String[], Float64[], Float64[])

	haskey(cfg, "reaction") || cfgerror("no [[reaction]] entries")
	rxns_raw = cfg["reaction"]
	rxns_raw isa AbstractVector || cfgerror("reactions must be written as [[reaction]] tables")
	reactions = ReactionSpec[]
	for (j, r) in enumerate(rxns_raw)
		check_keys(r, REACTION_KEYS, "[[reaction]] number $j")
		rname = string(get(r, "name", "R$j"))
		haskey(r, "reactant") || cfgerror("[[reaction]] $rname has no reactant")
		reactant = component_index(comps, string(r["reactant"]), "reaction $rname")
		for key in ("log10A", "E")
			haskey(r, key) || cfgerror("reaction $rname has no $key")
		end
		log10A = parse_slot!(fl, r["log10A"], "$rname.log10A")
		E = parse_slot!(fl, r["E"], "$rname.E"; lo = 0.0)
		order = parse_slot!(fl, get(r, "order", 1.0), "$rname.order"; lo = 0.0)
		prods_raw = get(r, "products", Dict{String, Any}())
		prods_raw isa AbstractDict || cfgerror("$rname.products must be a table, e.g. { char = 0.2 }")
		pnames = collect(keys(prods_raw))
		pidx = [component_index(comps, pn, "reaction $rname") for pn in pnames]
		products = Int[]
		yields = Slot[]
		for k in sortperm(pidx)
			p, pname = pidx[k], pnames[k]
			p == reactant && cfgerror("reaction $rname lists its reactant \"$pname\" as a product")
			push!(products, p)
			push!(yields, parse_slot!(fl, prods_raw[pname], "$rname.yield.$pname"; lo = 0.0, hi = 1.0))
		end
		fixed_sum = sum((y.value for y in yields if y.idx == 0); init = 0.0)
		fixed_sum <= 1.0 + 1e-12 || cfgerror("fixed yields of reaction $rname add up to $fixed_sum > 1")
		push!(reactions, ReactionSpec(rname, reactant, products, log10A, E, order, yields))
	end
	allunique(r.name for r in reactions) || cfgerror("reaction names must be unique")

	initial = [Slot(0, 0.0) for _ in comps]
	rest = 0
	init_raw = get(cfg, "initial", nothing)
	if init_raw === nothing
		initial[1] = Slot(0, 1.0)
	else
		init_raw isa AbstractDict || cfgerror("[initial] must be a table")
		for (cname, v) in init_raw
			i = component_index(comps, cname, "[initial]")
			if v == "rest"
				rest == 0 || cfgerror("only one [initial] entry can be \"rest\"")
				rest = i
			else
				initial[i] = parse_slot!(fl, v, "initial.$cname"; lo = 0.0, hi = 1.0)
			end
		end
		any_free = any(s -> s.idx != 0, initial)
		fixed_sum = sum((s.value for (i, s) in enumerate(initial) if s.idx == 0 && i != rest); init = 0.0)
		if rest == 0
			any_free && cfgerror("fitted [initial] fractions need one component set to \"rest\" so they add up to 1")
			abs(fixed_sum - 1.0) < 1e-6 || cfgerror("[initial] fractions add up to $fixed_sum, not 1 (or set one to \"rest\")")
		else
			fixed_sum <= 1.0 + 1e-12 || cfgerror("fixed [initial] fractions add up to $fixed_sum > 1")
		end
	end

	# A component that starts at zero and is never produced plays no part.
	for (i, c) in enumerate(comps)
		starts = i == rest || initial[i].idx != 0 || initial[i].value > 0
		made = any(r -> i in r.products, reactions)
		starts || made || @warn "component \"$c\" starts at zero and no reaction makes it"
	end

	return Scheme(comps, reactions, initial, rest, fl.names, fl.lower, fl.upper)
end

function describe(io::IO, s::Scheme)
	println(io, "Components: ", join(s.components, ", "))
	for r in s.reactions
		prods = isempty(r.products) ? "gas" :
				join(("$(slotstr(y)) $(s.components[p])" for (p, y) in zip(r.products, r.yields)), " + ") * " + gas"
		println(io, "  ", r.name, ": ", s.components[r.reactant], " -> ", prods)
	end
	println(io, "Free parameters: ", n_free(s))
end

slotstr(s::Slot) = s.idx == 0 ? @sprintf("%.4g", s.value) : "ν"
