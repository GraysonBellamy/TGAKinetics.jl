# TGAKinetics

Fit solid-phase reaction schemes to thermogravimetric (TGA) data measured in
nitrogen, at one or more heating rates at once. You describe the data files,
the reaction scheme and which parameters to fit in one text file; the tool
finds the Arrhenius parameters, solid yields (and, if you want, reaction
orders and initial composition) that best reproduce the measured mass and
mass-loss rate.

The fitting method is the one used in

> G. T. Bellamy, M. B. McKinnon, S. I. Stoliarov, Thermal decomposition models
> for wood: Comparing parallel and sequential reaction schemes for fire
> modeling applications, *Applications in Energy and Combustion Science* 26
> (2026) 100500. https://doi.org/10.1016/j.jaecs.2026.100500

The reaction form is the one FDS uses for solid-phase reactions: the fitted
A, E, order and yields map onto FDS's `A`, `E`, `N_S` and `NU_MATL` (convert
log10A to A, and kJ/mol to J/mol).

## Install

1. Install Julia 1.11 or newer. The simplest way is
   [juliaup](https://github.com/JuliaLang/juliaup):
   `curl -fsSL https://install.julialang.org | sh` (Linux/macOS) or
   `winget install julia -s msstore` (Windows).
2. Put this folder anywhere and, inside it, run once:

   ```
   julia --project=. -e "using Pkg; Pkg.instantiate()"
   ```

   This downloads the dependencies and compiles them (a few minutes).

You do not need to write any Julia. Everything is run from the command line
with a run file.

## Quick start

The folder `examples/synthetic` holds three made-up TGA runs (5, 10 and
20 K/min) generated from a known five-step scheme, so you can see what a
fit looks like before using your own data.

```
cd examples/synthetic
julia -t auto ../../bin/tgafit.jl check run.toml
julia -t auto ../../bin/tgafit.jl fit   run.toml
```

`-t auto` uses all CPU cores; `-t 4` uses four.

`check` reads the data and the scheme, reports what it found (reference mass,
peak mass-loss rate, number of compared points), writes the processed data to
`results/data_<name>.csv`, and estimates how long `fit` will take. Run it
first, and look at the processed data, before starting a long fit.

To plot the results (needs Python with numpy and matplotlib):

```
python ../../scripts/plot_results.py results/          # fit
python ../../scripts/plot_results.py results/ --data   # processed data from check
```

## Your data

One CSV file per TGA run, with a header row. Three columns are needed: time,
sample temperature and mass. Name them in the run file:

```toml
[data]
time_column = "Time (min)"
temperature_column = "Temperature (C)"
mass_column = "Mass (mg)"
time_unit = "min"          # "s" or "min"
temperature_unit = "C"     # "K" or "C"
```

If the header is not on the first line (instrument exports often start with
metadata), set `header_row` to its line number. Mass can be in any unit (it
is normalized). Use the measured sample
temperature, not the furnace set point; the heating program does not have to
be linear, and isothermal holds are fine.

What the tool does with each file:

1. **Mass-loss rate.** A cubic P-spline (500 segments, smoothing chosen by the
   V-curve) is fitted to mass against time, and its derivative is the
   mass-loss rate. If you already have a smoothed mass-loss rate, give it with
   `mlr_column` (as -dm/dt in the file's mass and time units).
2. **Reference mass.** With `reference = "dry"` (default), the reference is the
   smoothed mass at the point of lowest mass-loss rate inside `dry_window`
   (default 373-423 K), which removes adsorbed water. Data before that point
   are dropped and the model starts there. With `reference = "initial"` the
   first point is the reference; use this if you model moisture as a reaction.
3. **Temperature.** The measured temperature is smoothed with the same kind
   of spline, so thermocouple noise does not reach the model.
4. **Compared points.** Points are taken about every `grid_step` K (default 1),
   and only those inside `fit_range` (K) count in the objective.

All of these can be set in `[data]` for every file, or inside one
`[[dataset]]` to override them for that file.

```toml
[[dataset]]
name = "10Kmin"
file = "data/tga_10Kmin.csv"     # relative to the run file
weight = 0.5                      # weights are normalized to add up to 1
```

Replicates can be listed as separate datasets (for example each with weight
1/n of their heating rate's share) or averaged beforehand.

## Writing a reaction scheme

List the components, the initial composition and the reactions. Every value
can be **a number** (fixed) or **`[lower, upper]`** (fitted within those
bounds).

```toml
components = ["wood", "int", "char"]   # must come before the first [table]

[initial]            # mass fractions at the reference point, adding up to 1
wood = 1.0           # components not listed start at 0

[[reaction]]
name = "R1"
reactant = "wood"
products = { int = [0.5, 0.95] }   # solid yield; the rest goes to gas
log10A = [0.0, 18.0]               # log10 of A in 1/s
E = [40.0, 250.0]                  # kJ/mol
order = 1.0                        # optional, default 1

[[reaction]]
name = "R2"
reactant = "int"
products = { char = [0.1, 0.9] }
log10A = [0.0, 18.0]
E = [40.0, 250.0]
```

Each reaction consumes one reactant at the rate

  r = 10^log10A · exp(-E / (R T)) · (mass of reactant)^order

and turns it into the listed solid products (yield = mass of product per mass
of reactant consumed) plus gas (whatever the yields leave over). A reaction
with no `products` goes entirely to gas. A reaction can have several solid
products; if their fitted yields add up to more than 1, that candidate is
rejected.

What you can build this way:

- **Sequential schemes:** each reaction's product is the next reaction's
  reactant (`examples/synthetic/run.toml`).
- **Parallel (independent) components:** several components present at the
  start, each with its own reaction. Fit their initial fractions by giving
  ranges and setting one component to `"rest"` (`examples/synthetic/run_parallel.toml`):

  ```toml
  [initial]
  hemi = [0.1, 0.5]
  cell = [0.2, 0.6]
  lignin = "rest"     # 1 minus the others
  ```
- **Competing reactions:** two reactions with the same reactant.
- **Moisture:** a `water` component with a fixed initial fraction and a
  reaction with no products, together with `reference = "initial"`.
- **Inert residue:** a component with no reaction and a fixed initial fraction.
- **Known kinetics:** give numbers instead of ranges for any reaction you do
  not want to fit.

Typos in key names are reported as errors rather than ignored.

## Objective

The fit minimizes the composite normalized error (CNE) of Bellamy et al.
(2026), averaged over datasets with their weights:

  CNE = mlr_weight · MAE(MLR) / max(MLR measured) + mass_weight · MAE(mass) / max(mass measured)

where MAE is the mean absolute difference between model and measurement over
the compared points. 0 is a perfect fit. Defaults are 0.6 and 0.4:

```toml
[objective]
mlr_weight = 0.6
mass_weight = 0.4
```

As a guide, the CNE between replicate runs of the same wood at the same
heating rate was about 0.004-0.009 in the paper above; a fit cannot be
expected to do better than the scatter between your own replicates.

## Optimizer

Self-adaptive differential evolution (each population member keeps its own
control settings and mutation strategy), started from a Sobol sequence. The
population is scored in parallel on all threads.

```toml
[optimizer]
population = 50
generations = 1000
runs = 1          # independent repeats; the best one is kept
seed = 1          # same seed and settings give the same result
report_every = 10
ftol = 1e-8       # stop when the population's spread in CNE falls below this
xtol = 1e-8       # ... or its spread in the parameters
```

Fits with many parameters often settle in a local minimum, so one run is not
enough for a final answer. Use several `runs` (each starts from a different
part of the Sobol sequence with a different seed) and compare their results
in `runs.csv`. The paper used population 128, up to 100,000 generations,
`ftol = xtol = 1e-6`, and the best of 64 runs. Use `check` to see what a
setting will cost: one evaluation of three heating rates takes about 0.5 ms
here, so 128 × 20,000 is about 10 minutes on four cores per run.

ODE tolerances during the fit can be changed under `[solver]` (`abstol`,
default 1e-6; `reltol`, 1e-4; `maxiters`, 100000). Final curves and the
reported CNE are recomputed with 1e-9 / 1e-7.

## Output

In `output_dir` (default `results/`, next to the run file):

| File | Contents |
|---|---|
| `summary.txt` | CNE (total and per dataset) and the fitted parameters; parameters within 0.5 % of a bound are flagged |
| `parameters.csv` | fitted parameters with their bounds |
| `fitted.toml` | a complete run file with the fitted values filled in |
| `fit_<name>.csv` | time, temperature, measured and model mass and mass-loss rate, and every component's mass |
| `convergence.csv` | best CNE against generation |
| `runs.csv` | best CNE of each run (when `runs > 1`) |
| `best_so_far.csv` | written during a fit, removed at the end |

A parameter flagged as sitting on a bound usually means the bound should be
widened, or the scheme has more freedom than the data can pin down.

## Simulate

`simulate` runs a scheme with fixed values, for example the fitted one, on
data files or on set heating programs:

```
julia ../../bin/tgafit.jl simulate results/fitted.toml
```

To predict other heating rates, add datasets without a file:

```toml
[[dataset]]
name = "50Kmin"
heating_rate = 50.0    # K/min
T_start = 400.0        # K
T_end = 1100.0
```

Curves go to `sim_<name>.csv` in that run file's `output_dir`.

## Using it from Julia

```julia
using TGAKinetics
check("run.toml"); res = fit("run.toml"); simulate("results/fitted.toml")
```

## Tests

```
julia --project=. -e "using Pkg; Pkg.test()"
```

The tests compare the solver with exact solutions (first order, second order,
parallel components, linear heating), check the Jacobian against finite
differences, and run a short fit of the example end to end.

## Limits

- Nitrogen only: there is no oxygen dependence in the rate law.
- Mass only: heat flow (DSC) is not fitted.
- Reactions are single-reactant Arrhenius steps, as in FDS. Orders other than 1
  are supported but add stiffness near full conversion.
