using PowerDynamics
using OrderedCollections

const PLOT_TIMESERIES = false
const BASIN_PLOT = true

include("components.jl")

function findsym(pg)
    idx = findall(:ω .∈ values(pg.nodes) .|> symbolsof)
    return collect(keys(pg.nodes))[idx]
end

##
# Data Source: Kodsi, S. K. M., & Canizares, C. A. (2003). Modeling and simulation of IEEE 14-bus system with FACTS controllers. University of Waterloo, Canada, Tech. Rep.

# Assumptions:
# - select bus2 as the slack bus
# - replace generators with 3rd order model since we do not discuss higher-order SMs
# - loads will be interfaced with inverters

branches=OrderedDict(
    "branch1"=> PiModelLine(from= "bus1", to = "bus2",y=4.999131600798035-1im*15.263086523179553, y_shunt_km=0.0528/2, y_shunt_mk=0.0528/2),
    "branch2"=> PiModelLine(from= "bus1", to = "bus5",y=1.025897454970189-1im*4.234983682334831, y_shunt_km=0.0492/2, y_shunt_mk=0.0492/2),
    "branch3"=> PiModelLine(from= "bus2", to = "bus3",y=1.1350191923073958-1im*4.781863151757718, y_shunt_km=0.0438/2, y_shunt_mk=0.0438/2),
    "branch4"=> PiModelLine(from= "bus2", to = "bus4",y=1.686033150614943-1im*5.115838325872083, y_shunt_km=0.0340/2, y_shunt_mk=0.0340/2),
    "branch5"=> PiModelLine(from= "bus2", to = "bus5",y=1.7011396670944048-1im*5.193927397969713, y_shunt_km=0.0346/2, y_shunt_mk=0.0346/2),
    "branch6"=> PiModelLine(from= "bus3", to = "bus4",y=1.9859757099255606-1im*5.0688169775939205, y_shunt_km=0.0128/2, y_shunt_mk=0.0128/2),
    "branch7"=> StaticLine(from= "bus4", to = "bus5",Y=6.840980661495672-1im*21.578553981691588),
    "branch8"=> Transformer(from= "bus4", to = "bus7", y=0.0-1im*4.781943381790359, t_ratio=0.978),
    "branch9"=> Transformer(from= "bus4", to = "bus9", y=0.0-1im*1.7979790715236075, t_ratio=0.969),
    "branch10"=> Transformer(from= "bus5", to = "bus6", y=0.0-1im*3.967939052456154, t_ratio=0.932),
    "branch11"=> StaticLine(from= "bus6", to = "bus11",Y=1.9550285631772604-1im*4.0940743442404415),
    "branch12"=> StaticLine(from= "bus6", to = "bus12",Y=1.525967440450974-1im*3.1759639650294003),
    "branch13"=> StaticLine(from= "bus6", to = "bus13",Y=3.0989274038379877-1im*6.102755448193116),
    "branch14"=> StaticLine(from= "bus7", to = "bus8",Y=0.0-1im*5.676979846721544),
    "branch15"=> StaticLine(from= "bus7", to = "bus9",Y=0.0-1im*9.09008271975275),
    "branch16"=> StaticLine(from= "bus9", to = "bus10",Y=3.902049552447428-1im*10.365394127060915),
    "branch17"=> StaticLine(from= "bus9", to = "bus14",Y=1.4240054870199312-1im*3.0290504569306034),
    "branch18"=> StaticLine(from= "bus10", to = "bus11",Y=1.8808847537003996-1im*4.402943749460521),
    "branch19"=> StaticLine(from= "bus12", to = "bus13",Y=2.4890245868219187-1im*2.251974626172212),
    "branch20"=> StaticLine(from= "bus13", to = "bus14",Y=1.1369941578063267-1im*2.314963475105352));

# Schiffer models
τ_P=1
τ_Q=1
K_P=1
K_Q=1
# dVOC
η=0.001
α=1

buses_orig = OrderedDict(
    #"bus1"=> FourthOrderEq(T_d_dash=7.4, D=2, X_d=0.8979, X_q=0.646, Ω=50, X_d_dash=0.2995, T_q_dash=0.1, X_q_dash=0.646, P=2.32, H=5.148, E_f=1),
    "bus1"=> SchmietendorfOriginal(P_m=2.32, γ=2sqrt(50/(2*5.148)), α=7.4sqrt(50/(2*5.148)), E_f=1, X=0.8979),
    "bus2"=> SlackAlgebraic(U=1),
    #"bus3"=> FourthOrderEq(T_d_dash=6.1, D=2, X_d=1.05, X_q=0.98, Ω=50, X_d_dash=0.185, T_q_dash=0.4, X_q_dash=0.36, P=-0.942, H=6.54, E_f= 1),
    "bus3"=> SchmietendorfOriginal(P_m=-0.942, γ=2sqrt(50/(2*6.54)), α=6.1sqrt(50/(2*6.54)), E_f=1, X=1.05),
    #"bus4"=> VoltageDependentLoad(P=-0.478, Q=-0.0, U=1.0, A=0.0, B=0.0),
    "bus4" => SchifferOriginal(τ_P=τ_P, τ_Q=τ_Q, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.478, Q=0.0),
    #"bus5"=> VoltageDependentLoad(P=-0.076, Q=-0.016, U=1.0, A=0.0, B=0.0),
    "bus5" => SchifferOriginal(τ_P=τ_P, τ_Q=τ_Q, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.076, Q=-0.016),
    #"bus6"=> FourthOrderEq(T_d_dash=4.75, D=2, X_d=1.25, X_q=1.22, Ω=50, X_d_dash=0.232, T_q_dash=1.6, X_q_dash=0.715, P=-0.122, H=5.06, E_f= 1),
    "bus6"=> SchmietendorfOriginal(P_m=-0.122, γ=2sqrt(50/(2*5.06)), α=4.75sqrt(50/(2*5.06)), E_f=1, X=1.25),
    "bus7"=> PQAlgebraic(P=0.0, Q=0.0),
    #"bus8"=> FourthOrderEq(T_d_dash=4.75, D=2, X_d=1.25, X_q=1.22, Ω=50, X_d_dash=0.232, T_q_dash=1.6, X_q_dash=0.715, P=0.0, H=5.06, E_f= 1),
    "bus8"=> SchmietendorfOriginal(P_m=0., γ=2sqrt(50/(2*5.06)), α=4.75sqrt(50/(2*5.06)), E_f=1, X=1.25),
    #"bus9"=> VoltageDependentLoad(P=-0.295, Q=-0.166, U=1.0, A=0.0, B=0.0),
    "bus9" => SchifferOriginal(τ_P=τ_P, τ_Q=τ_Q, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.295, Q=-0.166),
    #"bus10"=> VoltageDependentLoad(P=-0.09, Q=-0.058, U=1.0, A=0.0, B=0.0),
    "bus10" => dVOC(Pset=-0.09, Qset=-0.058, Vset=1, Ωset=0, η=η, α=α, κ=π/2),
    #"bus11"=> VoltageDependentLoad(P=-0.035, Q=-0.018, U=1.0, A=0.0, B=0.0),
    "bus11"=> SchifferVIOriginal(;τ_P=τ_P, τ_Q=τ_Q, τ_V=0.005, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.035, Q=-0.018),
    #"bus12"=> VoltageDependentLoad(P=-0.061, Q=-0.016, U=1.0, A=0.0, B=0.0),
    "bus12" => SchifferOriginal(τ_P=τ_P, τ_Q=τ_Q, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.061, Q=-0.016),
    #"bus13"=> VoltageDependentLoad(P=-0.135, Q=-0.058, U=1.0, A=0.0, B=0.0),
    "bus13" => dVOC(Pset=-0.135, Qset=-0.058, Vset=1, Ωset=0, η=η, α=α, κ=π/2),
    #"bus13" => SchifferOriginal(τ_P=τ_P, τ_Q=τ_Q, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.135, Q=-0.058),
    #"bus14"=> VoltageDependentLoad(P=-0.149, Q=-0.05, U=1.0, A=0.0, B=0.0)
    "bus14"=> SchifferVIOriginal(;τ_P=τ_P, τ_Q=τ_Q, τ_V=0.005, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.149, Q=-0.05),
    );

buses_approx = OrderedDict(
    #"bus1"=> FourthOrderEq(T_d_dash=7.4, D=2, X_d=0.8979, X_q=0.646, Ω=50, X_d_dash=0.2995, T_q_dash=0.1, X_q_dash=0.646, P=2.32, H=5.148, E_f=1),
    "bus1"=> SchmietendorfApprox(P_m=2.32, γ=2sqrt(50/(2*5.148)), α=7.4sqrt(50/(2*5.148)), E_f=1, X=0.8979, E_c=1, Q_c=0),
    "bus2"=> SlackAlgebraic(U=1),
    #"bus3"=> FourthOrderEq(T_d_dash=6.1, D=2, X_d=1.05, X_q=0.98, Ω=50, X_d_dash=0.185, T_q_dash=0.4, X_q_dash=0.36, P=-0.942, H=6.54, E_f= 1),
    "bus3"=> SchmietendorfApprox(P_m=-0.942, γ=2sqrt(50/(2*6.54)), α=6.1sqrt(50/(2*6.54)), E_f=1, X=1.05, E_c=1, Q_c=0),
    #"bus4"=> VoltageDependentLoad(P=-0.478, Q=-0.0, U=1.0, A=0.0, B=0.0),
    "bus4" => SchifferApprox(τ_P=τ_P, τ_Q=τ_Q, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.478, Q=0.0),
    #"bus5"=> VoltageDependentLoad(P=-0.076, Q=-0.016, U=1.0, A=0.0, B=0.0),
    "bus5" => SchifferApprox(τ_P=τ_P, τ_Q=τ_Q, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.076, Q=-0.016),
    #"bus6"=> FourthOrderEq(T_d_dash=4.75, D=2, X_d=1.25, X_q=1.22, Ω=50, X_d_dash=0.232, T_q_dash=1.6, X_q_dash=0.715, P=-0.122, H=5.06, E_f= 1),
    "bus6"=> SchmietendorfApprox(P_m=-0.122, γ=2sqrt(50/(2*5.06)), α=4.75sqrt(50/(2*5.06)), E_f=1, X=1.25, E_c=1, Q_c=0),
    "bus7"=> PQAlgebraic(P=0.0, Q=0.0),
    #"bus8"=> FourthOrderEq(T_d_dash=4.75, D=2, X_d=1.25, X_q=1.22, Ω=50, X_d_dash=0.232, T_q_dash=1.6, X_q_dash=0.715, P=0.0, H=5.06, E_f= 1),
    "bus8"=> SchmietendorfApprox(P_m=0., γ=2sqrt(50/(2*5.06)), α=4.75sqrt(50/(2*5.06)), E_f=1, X=1.25, E_c=1, Q_c=0),
    #"bus9"=> VoltageDependentLoad(P=-0.295, Q=-0.166, U=1.0, A=0.0, B=0.0),
    "bus9" => SchifferApprox(τ_P=τ_P, τ_Q=τ_Q, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.295, Q=-0.166),
    #"bus10"=> VoltageDependentLoad(P=-0.09, Q=-0.058, U=1.0, A=0.0, B=0.0),
    "bus10" => dVOCapprox(Pset=-0.09, Qset=-0.058, Vset=1, Ωset=0, η=η, α=α, κ=π/2),
    #"bus11"=> VoltageDependentLoad(P=-0.035, Q=-0.018, U=1.0, A=0.0, B=0.0),
    "bus11"=> SchifferVIApprox(;τ_P=τ_P, τ_Q=τ_Q, τ_V=0.005, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.035, Q=-0.018),
    #"bus12"=> VoltageDependentLoad(P=-0.061, Q=-0.016, U=1.0, A=0.0, B=0.0),
    "bus12" => SchifferApprox(τ_P=τ_P, τ_Q=τ_Q, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.061, Q=-0.016),
    #"bus13"=> VoltageDependentLoad(P=-0.135, Q=-0.058, U=1.0, A=0.0, B=0.0),
    "bus13" => dVOCapprox(Pset=-0.135, Qset=-0.058, Vset=1, Ωset=0, η=η, α=α, κ=π/2),
    #"bus13" => SchifferApprox(τ_P=τ_P, τ_Q=τ_Q, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.135, Q=-0.058),
    #"bus14"=> VoltageDependentLoad(P=-0.149, Q=-0.05, U=1.0, A=0.0, B=0.0)
    "bus14"=> SchifferVIApprox(;τ_P=τ_P, τ_Q=τ_Q, τ_V=0.005, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.149, Q=-0.05),
    );

##

pg_orig = PowerGrid(buses_orig, branches)
op_orig = find_operationpoint(pg_orig; sol_method=:dynamic)

pg_approx = PowerGrid(buses_approx, branches)
op_approx = find_operationpoint(pg_approx; sol_method=:dynamic)

##

using Plots
using GraphPlot, Compose
using LaTeXStrings

plot_path = "$(@__DIR__)/../figures/"

ωidx = findsym(pg_orig)
##

println("finished setup")

##

if PLOT_TIMESERIES
    # @show op_orig[:, :φ];
    # @show op_approx[:, :φ];
    # @show op_orig[:, :v];
    # @show op_approx[:, :v];
    @show op_orig[ωidx, :ω];
    @show op_approx[ωidx, :ω];
    # @show op_orig[:, :p];
    # @show op_approx[:, :p];
    # @show op_orig[:, :q];
    # @show op_approx[:, :q];
    scatter(op_orig[:, :φ], op_approx[:, :φ], label="φ")
    scatter!(op_orig[:, :v], op_approx[:, :v], label="v")
    scatter!(op_orig[:, :p], op_approx[:, :p], label="p")
    scatter!(op_orig[:, :q], op_approx[:, :q], label="q")

    ##

    cmap = Dict(zip(
        [SchmietendorfOriginal, SlackAlgebraic, VSIMinimal, PQAlgebraic, dVOC, SchifferVIOriginal],
        distinguishable_colors(6, colorant"red")
    ))

    nc = map(x->cmap[x], pg_orig.nodes |> values .|> typeof)

    layout=(args...)->spring_layout(args...; C=20, INITTEMP=10)

    draw(SVG("$(plot_path)mixed_IEEE.svg"), gplot(pg_orig.graph, nodelabelc=RGB(1, 1, 1), nodefillc=nc, layout=layout,nodelabel=1:length(pg_orig.nodes)))


    ##

    tspan = (0., 50.)
    tarray = 0:0.01:last(tspan)
    node = "bus1"
    pp = PowerPerturbation(node=node, fault_power=1, tspan_fault=(10., 12.), var=:(P_m))

    sol_orig = simulate(pp, pg_orig, op_orig, tspan);
    sol_approx = simulate(pp, pg_approx, op_approx, tspan);

    ##

    plot(sol_orig, ["bus$i" for i in 1:5], :v, lc=:blue, legend=false, yguide="ρ", xguide="t")
    plot!(sol_orig, ["bus$i" for i in 6:14], :v, lc=:red, legend=false, yguide="ρ", xguide="t")

    ##

    plot(sol_orig, ["bus1", "bus3", "bus6", "bus8"], :v, lc=:blue, legend=false, yguide="ρ", xguide="t")
    plot!(sol_approx, ["bus1", "bus3", "bus6", "bus8"], :v, lc=:blue, ls=:dot)
    plot!(sol_orig, ["bus4", "bus5", "bus9", "bus12"], :v, lc=:green)
    plot!(sol_approx, ["bus4", "bus5", "bus9", "bus12"], :v, lc=:green, ls=:dot)
    plot!(sol_orig, ["bus10", "bus13"], :v, lc=:brown)
    plot!(sol_approx, ["bus10", "bus13"], :v, lc=:brown, ls=:dot)
    plot!(sol_orig, ["bus11", "bus14"], :v, lc=:orange)
    plot!(sol_approx, ["bus11", "bus14"], :v, lc=:orange, ls=:dot)
    vline!([pp.tspan_fault...], c=:black, alpha=0.8, label=false)

    savefig("$(plot_path)mixed_IEEE_trj_rho.png")

    ##

    plot(sol_orig, ["bus1", "bus3", "bus6", "bus8"], :φ, lc=:blue, legend=false, yguide="φ", xguide="t")
    plot!(sol_approx, ["bus1", "bus3", "bus6", "bus8"], :φ, lc=:blue, ls=:dot)
    plot!(sol_orig, ["bus4", "bus5", "bus9", "bus12"], :φ, lc=:green)
    plot!(sol_approx, ["bus4", "bus5", "bus9", "bus12"], :φ, lc=:green, ls=:dot)
    plot!(sol_orig, ["bus10", "bus13"], :φ, lc=:brown)
    plot!(sol_approx, ["bus10", "bus13"], :φ, lc=:brown, ls=:dot)
    plot!(sol_orig, ["bus11", "bus14"], :φ, lc=:orange)
    plot!(sol_approx, ["bus11", "bus14"], :φ, lc=:orange, ls=:dot)
    vline!([pp.tspan_fault...], c=:black, alpha=0.8, label=false)

    savefig("$(plot_path)mixed_IEEE_trj_phi.png")

    ##

    plot(sol_orig, ["bus1", "bus3", "bus6", "bus8"], :ω, lc=:blue, legend=false, yguide="ω", xguide="t")
    plot!(sol_approx, ["bus1", "bus3", "bus6", "bus8"], :ω, lc=:blue, ls=:dot)
    plot!(sol_orig, ["bus4", "bus5", "bus9", "bus12"], :ω, lc=:green)
    plot!(sol_approx, ["bus4", "bus5", "bus9", "bus12"], :ω, lc=:green, ls=:dot)
    plot!(sol_orig, ["bus11", "bus14"], :ω, lc=:orange)
    plot!(sol_approx, ["bus11", "bus14"], :ω, lc=:orange, ls=:dot)
    vline!([pp.tspan_fault...], c=:black, alpha=0.8, label=false)

    savefig("$(plot_path)mixed_IEEE_trj_omega.png")
end

##
println("begin basin plot section")
##

using ProbabilisticStability
using OrdinaryDiffEq
using QuasiMonteCarlo
using CSV
using Distances
using DataFrames

#

# git clone git@github.com:luap-pik/ProbabilisticStability.git
# ]dev path/to/repo

# We identify convergence given the frequency vanishes and ignore phase shifts.

##
function pg_dist_orig(x, y; d = Euclidean())
    sx = State(pg_orig, x)
    sy = State(pg_orig, y)
    return d([sx[:, :v]; sx[ωidx, :ω]], [sy[:, :v]; sy[ωidx, :ω]])
end

function pg_dist_approx(x, y; d = Euclidean())
    sx = State(pg_approx, x)
    sy = State(pg_approx, y)
    return d([sx[:, :v]; sx[ωidx, :ω]], [sy[:, :v]; sy[ωidx, :ω]])
end

mapping = Dict(:ω => L"$\omega \;[s^{-1}]$", :v => L"$\rho \;[pu]$", :φ => L"$\phi$")

function plotbasin(x, y, path, orig, approx, op_orig, xrange, yrange, model, node, vars, dimensions, perturbations)
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
            c=[:white, :yellow, :red, :orange],
            label=false,
            legend=false,
            xlims = extrema(xrange),
            ylims = extrema(yrange),
            xguide=x,
            yguide=y,
            size=(600, 600),
            tickfont=(24, "times"),
            guidefont=(26, "times"),
            legendfont=(18, "times"),
            framestyle=:box)
    scatter!([op_orig[node, first(vars)], ], [op_orig[node, last(vars)], ], markershape=:cross, markersize=10, c=:black, label="fix point")
    savefig(path)
end


# basin plot

function calc(OVERWRITE, model, node, vars, dimensions, perturbations)
    lb = [op_orig[node, first(vars)] - first(perturbations), op_orig[node, last(vars)] - last(perturbations)]
    ub = [op_orig[node, first(vars)] + first(perturbations), op_orig[node, last(vars)] + last(perturbations)]

    xrange = range(first(lb), stop=first(ub), length=Int(sqrt(sample_size)))
    yrange = range(last(lb), stop=last(ub), length=Int(sqrt(sample_size)))

    raw_ics = zeros(length(dimensions), sample_size)
    #raw_ics = QuasiMonteCarlo.sample(sample_size, lb, ub, LatticeRuleSample()) #GridSample([0.01,0.01])
    raw_ics[1, :] = [[x y] for x in xrange for y in yrange] .|> first #phase
    raw_ics[2, :] = [[x y] for x in xrange for y in yrange] .|> last #freq

    PD_ics = zeros(length(dimensions), sample_size)
    if length(dimensions) == 2
        PD_ics[1, :] = raw_ics[2, :] .* cos.(raw_ics[1, :])
        PD_ics[2, :] = raw_ics[2, :] .* sin.(raw_ics[1, :])
    else
        PD_ics[3, :] = raw_ics[2, :]
        PD_ics[1, :] = op_orig[node, :v] .* cos.(raw_ics[1, :])
        PD_ics[2, :] = op_orig[node, :v] .* sin.(raw_ics[1, :])
    end
    base_state = copy(op_orig.vec)
    base_state[dimensions] .= 0.
    ics = ICset(base_state, PD_ics, dimensions)

    if OVERWRITE || !isfile(plot_path * "$(model)_basin_orig.csv")
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

        CSV.write(plot_path * "$(model)_basin_orig.csv", results_orig_freq_i)
    end

    results_orig_freq_i = CSV.read(plot_path * "$(model)_basin_orig.csv", DataFrame)

    ##

    lb = [op_approx[node, first(vars)] - first(perturbations), op_approx[node, last(vars)] - last(perturbations)]
    ub = [op_approx[node, first(vars)] + first(perturbations), op_approx[node, last(vars)] + last(perturbations)]

    xrange = range(first(lb), stop=first(ub), length=Int(sqrt(sample_size)))
    yrange = range(last(lb), stop=last(ub), length=Int(sqrt(sample_size)))

    raw_ics = zeros(length(dimensions), sample_size)
    #raw_ics = QuasiMonteCarlo.sample(sample_size, lb, ub, LatticeRuleSample()) #GridSample([0.01,0.01])
    raw_ics[1, :] = [[x y] for x in xrange for y in yrange] .|> first
    raw_ics[2, :] = [[x y] for x in xrange for y in yrange] .|> last

    PD_ics = zeros(length(dimensions), sample_size)
    if length(dimensions) == 2
        PD_ics[1, :] = raw_ics[2, :] .* cos.(raw_ics[1, :])
        PD_ics[2, :] = raw_ics[2, :] .* sin.(raw_ics[1, :])
    else
        PD_ics[3, :] = raw_ics[2, :]
        PD_ics[1, :] = op_approx[node, :v] .* cos.(raw_ics[1, :])
        PD_ics[2, :] = op_approx[node, :v] .* sin.(raw_ics[1, :])
    end
    base_state = copy(op_approx.vec)
    base_state[dimensions] .= 0.
    ics = ICset(base_state, PD_ics, dimensions)

    if OVERWRITE || !isfile(plot_path * "$(model)_basin_approx.csv")
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

        CSV.write(plot_path * "$(model)_basin_approx.csv", results_approx_freq_i)
    end

    results_approx_freq_i = CSV.read(plot_path * "$(model)_basin_approx.csv", DataFrame)

    ##

    plotbasin(
        mapping[first(vars)], mapping[last(vars)],
        plot_path * "$(model)_basin.png",
        results_orig_freq_i, results_approx_freq_i, op_orig, xrange, yrange,
        model, node, vars, dimensions, perturbations
        )
end

##
println("finished function defs")
##

dim = map(dimension, values(pg_orig.nodes))
last_idx =   dim |> cumsum
first_idx = 1 .+ last_idx .- dim;
@show zip.(first_idx, last_idx)[13];

##
setups = [
    ("bus1", [:φ, :v], 1:2, (π, 0.8)),
    ("bus4", [:φ, :ω], 9:11, (π, 10)),
    ("bus9", [:φ, :v], 23:24, (π, 0.8)),
    ("bus6", [:v, :ω], 15:17, (0.8, 10)),
    ("bus8", [:φ, :ω], 20:22, (π, 10)),
    ("bus12", [:φ, :v], 32:33, (π, 0.8)),
    ("bus13", [:φ, :v], 35:36, (π, 0.8))
]

sample_size = 16900 # square number!!
#vars = [:φ, :v]
#dimensions = 26:27 #map(dimension, values(pg_orig.nodes)) |> cumsum
#perturbations = (π, 0.8) # (π, 10) # (π, 0.8)
OVERWRITE = false
#node = "bus10"

if BASIN_PLOT
    for setup in setups
        node, vars, dimensions, perturbations = setup
        model = "mixed_IEEE" * "_" * node * "_" * string(first(vars)) * "_" * string(last(vars))
        calc(OVERWRITE, model, node, vars, dimensions, perturbations)
    end
end
