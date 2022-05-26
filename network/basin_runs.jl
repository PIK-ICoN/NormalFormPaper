using ProbabilisticStability
using OrdinaryDiffEq
using QuasiMonteCarlo
using CSV
using DataFrames


##

function plotbasin(x, y, path, orig, approx, op_orig, xrange, yrange)
    @assert size(orig, 1) == size(approx, 1)
    sample_size = size(orig, 1)
    z1 = map(x->x ? 1 : 0, orig.within_threshold .& (orig.retcode .!== "Unstable"))
    z2 = map(x->x ? 2 : 0, approx.within_threshold .& (approx.retcode .!== "Unstable"))
    heatmap(
            xrange,
            yrange,
            reshape(z1 .+ z2, Int(sqrt(sample_size)), :),
            levels=4,
            colorbar=false,
            grid=false,
            c=cgrad([:white, :yellow, :red, :orange]),
            label=false,
            legend=false,
            xlims = extrema(xrange),
            ylims = extrema(yrange),
            size=(600, 600),
            xguide=x,
            yguide=y,
            fontfamily="Computer Modern",
            tickfontsize=24,
            guidefontsize=24,
            legendfontsize=18,
            framestyle=:box)
    scatter!([op_orig[node, first(vars)], ], [op_orig[node, last(vars)], ], markershape=:cross, markersize=10, c=:black, label="fix point") |> display
    savefig(path)
end

## basin φ, u

node = 1
vars = [:φ, :v]
dimensions = 1:2

lb = [op_orig[node, first(vars)] - π, op_orig[node, last(vars)] - 1]
ub = [op_orig[node, first(vars)] + π, op_orig[node, last(vars)] + 1]
xrange = range(first(lb), stop=first(ub), length=Int(sqrt(sample_size)))
yrange = range(last(lb), stop=last(ub), length=Int(sqrt(sample_size)))

raw_ics = zeros(length(dimensions), sample_size)
raw_ics[1, :] = [[x y] for x in xrange for y in yrange] .|> first #phase
raw_ics[2, :] = [[x y] for x in xrange for y in yrange] .|> last #amplitude

PD_ics = zeros(length(dimensions), sample_size)
PD_ics[1, :] = raw_ics[2, :] .* cos.(raw_ics[1, :])
PD_ics[2, :] = raw_ics[2, :] .* sin.(raw_ics[1, :])
base_state = copy(op_orig.vec)
base_state[dimensions] .= 0.
ics = ICset(base_state, PD_ics, dimensions)

if !isfile(plot_path * "$(model)_basin_rho_phi_orig.csv")
    @time (μ, μerr, total), results_orig = basin_stability_fixpoint(
        pg_orig,
        op_orig.vec,
        ics;
        distance = pg_dist_orig,
        threshold = 1E-4,
        parallel_alg = EnsembleThreads(),
        solver = Rodas5(),
        verbose = true,
        return_df = true,
        )

    CSV.write(plot_path * "$(model)_basin_rho_phi_orig.csv", results_orig)
end

results_orig = CSV.read(plot_path * "$(model)_basin_rho_phi_orig.csv", DataFrame)
##

node = 1
vars = [:φ, :v]
dimensions = 1:2

lb = [op_approx[node, first(vars)] - π, op_approx[node, last(vars)] - 1]
ub = [op_approx[node, first(vars)] + π, op_approx[node, last(vars)] + 1]
xrange = range(first(lb), stop=first(ub), length=Int(sqrt(sample_size)))
yrange = range(last(lb), stop=last(ub), length=Int(sqrt(sample_size)))

raw_ics = zeros(length(dimensions), sample_size)
raw_ics[1, :] = [[x y] for x in xrange for y in yrange] .|> first
raw_ics[2, :] = [[x y] for x in xrange for y in yrange] .|> last

PD_ics = zeros(length(dimensions), sample_size)
PD_ics[1, :] = raw_ics[2, :] .* cos.(raw_ics[1, :])

PD_ics[2, :] = raw_ics[2, :] .* sin.(raw_ics[1, :])
base_state = copy(op_approx.vec)
base_state[dimensions] .= 0.
ics = ICset(base_state, PD_ics, dimensions)

if !isfile(plot_path * "$(model)_basin_rho_phi_approx.csv")
    @time (μ, μerr, total), results_approx = basin_stability_fixpoint(
        pg_approx,
        op_approx.vec,
        ics;
        distance = pg_dist_approx,
        threshold = 1E-4,
        parallel_alg = EnsembleThreads(),
        solver = Rodas5(),
        verbose = true,
        return_df = true,
        )

    CSV.write(plot_path * "$(model)_basin_rho_phi_approx.csv", results_approx)
end

results_approx = CSV.read(plot_path * "$(model)_basin_rho_phi_approx.csv", DataFrame)

##

# basin_plot(op_orig, results_orig, raw_ics, lb, ub; node=1, syms=vars, ms=3, label=false, xguide=L"$\varphi$", yguide=L"$\rho$")
# scatter!([op_orig[node, first(vars)], ], [op_orig[node, last(vars)], ], markershape=:star, markersize=10, c=:black, label="fix point")
# basin_plot(op_approx, results_approx, raw_ics, lb, ub; node=1, syms=vars, ms=3, label=false, xguide=L"$\varphi$", yguide=L"$\rho$")
# scatter!([op_approx[node, first(vars)], ], [op_approx[node, last(vars)], ], markershape=:star, markersize=10, c=:black, label="fix point")

plotbasin(
    L"$\phi \;(rad)$", L"$\rho \;(pu)$",
    plot_path * "$(model)_basin_rho_phi.png",
    results_orig, results_approx, op_orig, xrange, yrange
    )

## basin u, ω

node = 1
vars = [:ω, :v]
dimensions = [1, 2, 3]

lb = [op_orig[node, first(vars)] - 10, op_orig[node, last(vars)] - 1]
ub = [op_orig[node, first(vars)] + 10, op_orig[node, last(vars)] + 1]

xrange = range(first(lb), stop=first(ub), length=Int(sqrt(sample_size)))
yrange = range(last(lb), stop=last(ub), length=Int(sqrt(sample_size)))

raw_ics = zeros(length(dimensions), sample_size)
#raw_ics = QuasiMonteCarlo.sample(sample_size, lb, ub, LatticeRuleSample()) #GridSample([0.01,0.01])
raw_ics[1, :] = [[x y] for x in xrange for y in yrange] .|> first # freq
raw_ics[2, :] = [[x y] for x in xrange for y in yrange] .|> last # amplitude

PD_ics = zeros(length(dimensions), sample_size)
PD_ics[3, :] = raw_ics[1, :]
PD_ics[1, :] = raw_ics[2, :] .* cos.(op_orig[node, :φ])
PD_ics[2, :] = raw_ics[2, :] .* sin.(op_orig[node, :φ])
base_state = copy(op_orig.vec)
base_state[dimensions] .= 0.
ics = ICset(base_state, PD_ics, dimensions)

if !isfile(plot_path * "$(model)_basin_rho_omega_orig.csv")
    @time (μ, μerr, total), results_orig_freq = basin_stability_fixpoint(
        pg_orig,
        op_orig.vec,
        ics;
        distance = pg_dist_orig,
        threshold = 1E-4,
        parallel_alg = EnsembleThreads(),
        solver = Rodas5(),
        verbose = true,
        return_df = true,
        )

    CSV.write(plot_path * "$(model)_basin_rho_omega_orig.csv", results_orig_freq)
end

results_orig_freq = CSV.read(plot_path * "$(model)_basin_rho_omega_orig.csv", DataFrame)


##

node = 1
vars = [:ω, :v]
dimensions = [1, 2, 3]

lb = [op_approx[node, first(vars)] - 10, op_approx[node, last(vars)] - 1]
ub = [op_approx[node, first(vars)] + 10, op_approx[node, last(vars)] + 1]

xrange = range(first(lb), stop=first(ub), length=Int(sqrt(sample_size)))
yrange = range(last(lb), stop=last(ub), length=Int(sqrt(sample_size)))

raw_ics = zeros(length(dimensions), sample_size)
#raw_ics = QuasiMonteCarlo.sample(sample_size, lb, ub, LatticeRuleSample()) #GridSample([0.01,0.01])
raw_ics[1, :] = [[x y] for x in xrange for y in yrange] .|> first
raw_ics[2, :] = [[x y] for x in xrange for y in yrange] .|> last

PD_ics = zeros(length(dimensions), sample_size)
PD_ics[3, :] = raw_ics[1, :]
PD_ics[1, :] = raw_ics[2, :] .* cos.(op_approx[node, :φ])
PD_ics[2, :] = raw_ics[2, :] .* sin.(op_approx[node, :φ])
base_state = copy(op_approx.vec)
base_state[dimensions] .= 0.
ics = ICset(base_state, PD_ics, dimensions)

if !isfile(plot_path * "$(model)_basin_rho_omega_approx.csv")
    @time (μ, μerr, total), results_approx_freq = basin_stability_fixpoint(
        pg_approx,
        op_approx.vec,
        ics;
        distance = pg_dist_approx,
        threshold = 1E-4,
        parallel_alg = EnsembleThreads(),
        solver = Rodas5(),
        verbose = true,
        return_df = true,
        )

    CSV.write(plot_path * "$(model)_basin_rho_omega_approx.csv", results_approx_freq)
end

results_approx_freq = CSV.read(plot_path * "$(model)_basin_rho_omega_approx.csv", DataFrame)

##

# basin_plot(op_orig, results_orig_freq, raw_ics, lb, ub; node=1, syms=vars, ms=3, label=false, xguide=L"$\omega$", yguide=L"$\rho$")
# scatter!([op_orig[node, first(vars)], ], [op_orig[node, last(vars)], ], markershape=:star, markersize=10, c=:black, label="fix point")
# basin_plot(op_approx, results_approx_freq, raw_ics, lb, ub; node=1, syms=vars, ms=3, label=false, xguide=L"$\omega$", yguide=L"$\rho$")
# scatter!([op_approx[node, first(vars)], ], [op_approx[node, last(vars)], ], markershape=:star, markersize=10, c=:black, label="fix point")

plotbasin(
    L"$\omega \;(rad/s)$", L"$\rho \;(pu)$",
    plot_path * "$(model)_basin_rho_omega.png",
    results_orig_freq, results_approx_freq, op_orig, xrange, yrange
    )

## basin ω, φ

node = 1
vars = [:φ, :ω]
dimensions = [1, 2, 3]

lb = [op_orig[node, first(vars)] - π, op_orig[node, last(vars)] - 10]
ub = [op_orig[node, first(vars)] + π, op_orig[node, last(vars)] + 10]

xrange = range(first(lb), stop=first(ub), length=Int(sqrt(sample_size)))
yrange = range(last(lb), stop=last(ub), length=Int(sqrt(sample_size)))

raw_ics = zeros(length(dimensions), sample_size)
#raw_ics = QuasiMonteCarlo.sample(sample_size, lb, ub, LatticeRuleSample()) #GridSample([0.01,0.01])
raw_ics[1, :] = [[x y] for x in xrange for y in yrange] .|> first #phase
raw_ics[2, :] = [[x y] for x in xrange for y in yrange] .|> last #freq

PD_ics = zeros(length(dimensions), sample_size)
PD_ics[3, :] = raw_ics[2, :]
PD_ics[1, :] = op_orig[node, :v] .* cos.(raw_ics[1, :])
PD_ics[2, :] = op_orig[node, :v] .* sin.(raw_ics[1, :])
base_state = copy(op_orig.vec)
base_state[dimensions] .= 0.
ics = ICset(base_state, PD_ics, dimensions)

if !isfile(plot_path * "$(model)_basin_omega_phi_orig.csv")
    @time (μ, μerr, total), results_orig_freq_i = basin_stability_fixpoint(
        pg_orig,
        op_orig.vec,
        ics;
        distance = pg_dist_orig,
        threshold = 1E-4,
        parallel_alg = EnsembleThreads(),
        solver = Rodas5(),
        verbose = true,
        return_df = true,
        )

    CSV.write(plot_path * "$(model)_basin_omega_phi_orig.csv", results_orig_freq_i)
end

results_orig_freq_i = CSV.read(plot_path * "$(model)_basin_omega_phi_orig.csv", DataFrame)

##

node = 1
vars = [:φ, :ω]
dimensions = [1, 2, 3]

lb = [op_approx[node, first(vars)] - π, op_approx[node, last(vars)] - 10]
ub = [op_approx[node, first(vars)] + π, op_approx[node, last(vars)] + 10]

xrange = range(first(lb), stop=first(ub), length=Int(sqrt(sample_size)))
yrange = range(last(lb), stop=last(ub), length=Int(sqrt(sample_size)))

raw_ics = zeros(length(dimensions), sample_size)
#raw_ics = QuasiMonteCarlo.sample(sample_size, lb, ub, LatticeRuleSample()) #GridSample([0.01,0.01])
raw_ics[1, :] = [[x y] for x in xrange for y in yrange] .|> first
raw_ics[2, :] = [[x y] for x in xrange for y in yrange] .|> last

PD_ics = zeros(length(dimensions), sample_size)
PD_ics[3, :] = raw_ics[2, :]
PD_ics[1, :] = op_approx[node, :v] .* cos.(raw_ics[1, :])
PD_ics[2, :] = op_approx[node, :v] .* sin.(raw_ics[1, :])
base_state = copy(op_approx.vec)
base_state[dimensions] .= 0.
ics = ICset(base_state, PD_ics, dimensions)

if !isfile(plot_path * "$(model)_basin_omega_phi_approx.csv")
    @time (μ, μerr, total), results_approx_freq_i = basin_stability_fixpoint(
        pg_approx,
        op_approx.vec,
        ics;
        distance = pg_dist_approx,
        threshold = 1E-4,
        parallel_alg = EnsembleThreads(),
        solver = Rodas5(),
        verbose = true,
        return_df = true,
        )

    CSV.write(plot_path * "$(model)_basin_omega_phi_approx.csv", results_approx_freq_i)
end

results_approx_freq_i = CSV.read(plot_path * "$(model)_basin_omega_phi_approx.csv", DataFrame)

##

# basin_plot(op_orig, results_orig_freq_i, raw_ics, lb, ub; node=1, syms=vars, ms=3, label=false, xguide=L"$\varphi$", yguide=L"$\omega$")
# scatter!([op_orig[node, first(vars)], ], [op_orig[node, last(vars)], ], markershape=:star, markersize=10, c=:black, label="fix point")
# basin_plot(op_approx, results_approx_freq_i, raw_ics, lb, ub; node=1, syms=vars, ms=3, label=false, xguide=L"$\varphi$", yguide=L"$\omega$")
# scatter!([op_approx[node, first(vars)], ], [op_approx[node, last(vars)], ], markershape=:star, markersize=10, c=:black, label="fix point")

plotbasin(
    L"$\phi \;(rad)$", L"$\omega \;(rad/s)$",
    plot_path * "$(model)_basin_omega_phi.png",
    results_orig_freq_i, results_approx_freq_i, op_orig, xrange, yrange
    )
