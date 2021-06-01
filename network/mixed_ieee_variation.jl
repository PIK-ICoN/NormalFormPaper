using PowerDynamics
using OrderedCollections
using Distances

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

##
# Schiffer models
function original_model(;scaling=1, τ_P=1, τ_Q=1, K_P=1, K_Q=1, η=0.001, α=1)
    buses_orig = OrderedDict(
        #"bus1"=> FourthOrderEq(T_d_dash=7.4, D=2, X_d=0.8979, X_q=0.646, Ω=50, X_d_dash=0.2995, T_q_dash=0.1, X_q_dash=0.646, P=2.32, H=5.148, E_f=1),
        "bus1"=> SchmietendorfOriginal(P_m=2.32, γ=2sqrt(50/(2*5.148)), α=7.4sqrt(50/(2*5.148)), E_f=1, X=0.8979),
        "bus2"=> SlackAlgebraic(U=1),
        #"bus3"=> FourthOrderEq(T_d_dash=6.1, D=2, X_d=1.05, X_q=0.98, Ω=50, X_d_dash=0.185, T_q_dash=0.4, X_q_dash=0.36, P=-0.942, H=6.54, E_f= 1),
        "bus3"=> SchmietendorfOriginal(P_m=-0.942, γ=2sqrt(50/(2*6.54)), α=6.1sqrt(50/(2*6.54)), E_f=1, X=1.05),
        #"bus4"=> VoltageDependentLoad(P=-0.478, Q=-0.0, U=1.0, A=0.0, B=0.0),
        "bus4" => SchifferOriginal(τ_P=τ_P, τ_Q=τ_Q, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.478*scaling, Q=0.0*scaling),
        #"bus5"=> VoltageDependentLoad(P=-0.076, Q=-0.016, U=1.0, A=0.0, B=0.0),
        "bus5" => SchifferOriginal(τ_P=τ_P, τ_Q=τ_Q, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.076*scaling, Q=-0.016*scaling),
        #"bus6"=> FourthOrderEq(T_d_dash=4.75, D=2, X_d=1.25, X_q=1.22, Ω=50, X_d_dash=0.232, T_q_dash=1.6, X_q_dash=0.715, P=-0.122, H=5.06, E_f= 1),
        "bus6"=> SchmietendorfOriginal(P_m=-0.122, γ=2sqrt(50/(2*5.06)), α=4.75sqrt(50/(2*5.06)), E_f=1, X=1.25),
        "bus7"=> PQAlgebraic(P=0.0, Q=0.0),
        #"bus8"=> FourthOrderEq(T_d_dash=4.75, D=2, X_d=1.25, X_q=1.22, Ω=50, X_d_dash=0.232, T_q_dash=1.6, X_q_dash=0.715, P=0.0, H=5.06, E_f= 1),
        "bus8"=> SchmietendorfOriginal(P_m=0., γ=2sqrt(50/(2*5.06)), α=4.75sqrt(50/(2*5.06)), E_f=1, X=1.25),
        #"bus9"=> VoltageDependentLoad(P=-0.295, Q=-0.166, U=1.0, A=0.0, B=0.0),
        "bus9" => SchifferOriginal(τ_P=τ_P, τ_Q=τ_Q, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.295*scaling, Q=-0.166*scaling),
        #"bus10"=> VoltageDependentLoad(P=-0.09, Q=-0.058, U=1.0, A=0.0, B=0.0),
        "bus10" => dVOC(Pset=-0.09*scaling, Qset=-0.058*scaling, Vset=1, Ωset=0, η=η, α=α, κ=π/2),
        #"bus11"=> VoltageDependentLoad(P=-0.035, Q=-0.018, U=1.0, A=0.0, B=0.0),
        "bus11"=> SchifferVIOriginal(;τ_P=τ_P, τ_Q=τ_Q, τ_V=0.005, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.035*scaling, Q=-0.018*scaling),
        #"bus12"=> VoltageDependentLoad(P=-0.061, Q=-0.016, U=1.0, A=0.0, B=0.0),
        "bus12" => SchifferOriginal(τ_P=τ_P, τ_Q=τ_Q, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.061*scaling, Q=-0.016*scaling),
        #"bus13"=> VoltageDependentLoad(P=-0.135, Q=-0.058, U=1.0, A=0.0, B=0.0),
        "bus13" => dVOC(Pset=-0.135*scaling, Qset=-0.058*scaling, Vset=1, Ωset=0, η=η, α=α, κ=π/2),
        #"bus14"=> VoltageDependentLoad(P=-0.149, Q=-0.05, U=1.0, A=0.0, B=0.0)
        "bus14"=> SchifferVIOriginal(;τ_P=τ_P, τ_Q=τ_Q, τ_V=0.005, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.149*scaling, Q=-0.05*scaling),
        );

        pg = PowerGrid(buses_orig, branches)
        ωidx = findsym(pg)

        function pg_dist(x, y; d = Euclidean())
            sx = State(pg, x)
            sy = State(pg, y)
            return d([sx[:, :v]; sx[ωidx, :ω]], [sy[:, :v]; sy[ωidx, :ω]])
        end

    return pg, find_operationpoint(pg; sol_method=:dynamic), pg_dist
end

function approximate_model(;scaling=1, τ_P=1, τ_Q=1, K_P=1, K_Q=1, η=0.001, α=1)
    buses_approx = OrderedDict(
        #"bus1"=> FourthOrderEq(T_d_dash=7.4, D=2, X_d=0.8979, X_q=0.646, Ω=50, X_d_dash=0.2995, T_q_dash=0.1, X_q_dash=0.646, P=2.32, H=5.148, E_f=1),
        "bus1"=> SchmietendorfApprox(P_m=2.32, γ=2sqrt(50/(2*5.148)), α=7.4sqrt(50/(2*5.148)), E_f=1, X=0.8979, E_c=1, Q_c=0),
        "bus2"=> SlackAlgebraic(U=1),
        #"bus3"=> FourthOrderEq(T_d_dash=6.1, D=2, X_d=1.05, X_q=0.98, Ω=50, X_d_dash=0.185, T_q_dash=0.4, X_q_dash=0.36, P=-0.942, H=6.54, E_f= 1),
        "bus3"=> SchmietendorfApprox(P_m=-0.942, γ=2sqrt(50/(2*6.54)), α=6.1sqrt(50/(2*6.54)), E_f=1, X=1.05, E_c=1, Q_c=0),
        #"bus4"=> VoltageDependentLoad(P=-0.478, Q=-0.0, U=1.0, A=0.0, B=0.0),
        "bus4" => SchifferApprox(τ_P=τ_P, τ_Q=τ_Q, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.478*scaling, Q=0.0*scaling),
        #"bus5"=> VoltageDependentLoad(P=-0.076, Q=-0.016, U=1.0, A=0.0, B=0.0),
        "bus5" => SchifferApprox(τ_P=τ_P, τ_Q=τ_Q, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.076*scaling, Q=-0.016*scaling),
        #"bus6"=> FourthOrderEq(T_d_dash=4.75, D=2, X_d=1.25, X_q=1.22, Ω=50, X_d_dash=0.232, T_q_dash=1.6, X_q_dash=0.715, P=-0.122, H=5.06, E_f= 1),
        "bus6"=> SchmietendorfApprox(P_m=-0.122, γ=2sqrt(50/(2*5.06)), α=4.75sqrt(50/(2*5.06)), E_f=1, X=1.25, E_c=1, Q_c=0),
        "bus7"=> PQAlgebraic(P=0.0, Q=0.0),
        #"bus8"=> FourthOrderEq(T_d_dash=4.75, D=2, X_d=1.25, X_q=1.22, Ω=50, X_d_dash=0.232, T_q_dash=1.6, X_q_dash=0.715, P=0.0, H=5.06, E_f= 1),
        "bus8"=> SchmietendorfApprox(P_m=0., γ=2sqrt(50/(2*5.06)), α=4.75sqrt(50/(2*5.06)), E_f=1, X=1.25, E_c=1, Q_c=0),
        #"bus9"=> VoltageDependentLoad(P=-0.295, Q=-0.166, U=1.0, A=0.0, B=0.0),
        "bus9" => SchifferApprox(τ_P=τ_P, τ_Q=τ_Q, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.295*scaling, Q=-0.166*scaling),
        #"bus10"=> VoltageDependentLoad(P=-0.09, Q=-0.058, U=1.0, A=0.0, B=0.0),
        "bus10" => dVOCapprox(Pset=-0.09*scaling, Qset=-0.058*scaling, Vset=1, Ωset=0, η=η, α=α, κ=π/2),
        #"bus11"=> VoltageDependentLoad(P=-0.035, Q=-0.018, U=1.0, A=0.0, B=0.0),
        "bus11"=> SchifferVIApprox(;τ_P=τ_P, τ_Q=τ_Q, τ_V=0.005, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.035*scaling, Q=-0.018*scaling),
        #"bus12"=> VoltageDependentLoad(P=-0.061, Q=-0.016, U=1.0, A=0.0, B=0.0),
        "bus12" => SchifferApprox(τ_P=τ_P, τ_Q=τ_Q, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.061*scaling, Q=-0.016*scaling),
        #"bus13"=> VoltageDependentLoad(P=-0.135, Q=-0.058, U=1.0, A=0.0, B=0.0),
        "bus13" => dVOCapprox(Pset=-0.135*scaling, Qset=-0.058*scaling, Vset=1, Ωset=0, η=η, α=α, κ=π/2),
        #"bus14"=> VoltageDependentLoad(P=-0.149, Q=-0.05, U=1.0, A=0.0, B=0.0)
        "bus14"=> SchifferVIApprox(;τ_P=τ_P, τ_Q=τ_Q, τ_V=0.005, K_P=K_P, K_Q=K_Q, V_r=1, P=-0.149*scaling, Q=-0.05*scaling),
        );

    pg = PowerGrid(buses_approx, branches)
    ωidx = findsym(pg)

    function pg_dist(x, y; d = Euclidean())
        sx = State(pg, x)
        sy = State(pg, y)
        return d([sx[:, :v]; sx[ωidx, :ω]], [sy[:, :v]; sy[ωidx, :ω]])
    end
    return pg, find_operationpoint(pg; sol_method=:dynamic), pg_dist
end

##

pg_orig, op_orig, _ = original_model()
pg_approx, op_approx, _ = approximate_model()

##

using Plots
using GraphPlot, Compose
using LaTeXStrings

plot_path = "$(@__DIR__)/../figures/"


##

function stability_interval(sf)
    pg_orig, op_orig, _ = original_model(; scaling=sf)
    pg_approx, op_approx, _ = approximate_model(; scaling=sf)
    ωidx = findsym(pg_orig)
    return minimum(op_orig[:, :v]), minimum(op_approx[:, :v]), maximum(op_orig[ωidx,:ω]), maximum(op_approx[ωidx,:ω])
end

##

@show stability_interval(2);

#################################################################################################################################
## variation of power demand

# git clone git@github.com:luap-pik/ProbabilisticStability.git
# ]dev path/to/repo
using ProbabilisticStability
using Measurements
using OrdinaryDiffEq
using Distances
using QuasiMonteCarlo
using CSV, DataFrames


############################################################
# very long expression starts here
##

function dimensional_basin_analyses(pg, op, pg_dist)
    # perturb all dimension

    lb = [op[node, :φ] - perturbations[:φ], op[node, :v] - perturbations[:v], op[node, :ω] - perturbations[:ω]]
    ub = [op[node, :φ] + perturbations[:φ], op[node, :v] + perturbations[:v], op[node, :ω] + perturbations[:ω]]

    raw_ics = QuasiMonteCarlo.sample(sample_size, lb, ub, SobolSample())

    PD_ics = zeros(length(dimensions), sample_size)
    PD_ics[3, :] = raw_ics[3, :]
    PD_ics[1, :] = raw_ics[2, :] .* cos.(raw_ics[1, :])
    PD_ics[2, :] = raw_ics[2, :] .* sin.(raw_ics[1, :])
    base_state = copy(op.vec)
    base_state[dimensions] .= 0.
    ics = ICset(base_state, PD_ics, dimensions)


    @time (μ, μerr, total)= basin_stability_fixpoint(
        pg,
        op.vec,
        ics;
        distance = pg_dist,
        threshold = 1E-4,
        parallel_alg = EnsembleThreads(),
        solver = Rodas5(),
        verbose = true,
        return_df = false,
        )

    combined = μ ± μerr

    # # perturb frequency

    # lb = op[node, :ω] - perturbations[:ω]
    # ub = op[node, :ω] + perturbations[:ω]

    # raw_ics = QuasiMonteCarlo.sample(sample_size, lb, ub, SobolSample())
    # base_state = copy(op.vec)
    # base_state[3] = 0.
    # ics = ICset(base_state, raw_ics', 3)

    # @time (μ, μerr, total)= basin_stability_fixpoint(
    #     pg,
    #     op.vec,
    #     ics;
    #     distance = pg_dist,
    #     threshold = 1E-4,
    #     parallel_alg = EnsembleThreads(),
    #     solver = Rodas5(),
    #     verbose = true,
    #     return_df = false,
    #     )

    # freq = μ ± μerr

    # # perturb amplitude

    # lb = op[node, :v] - perturbations[:v]
    # ub = op[node, :v] + perturbations[:v]

    # raw_ics = QuasiMonteCarlo.sample(sample_size, lb, ub, SobolSample())

    # PD_ics = zeros(2, sample_size)
    # PD_ics[1, :] = raw_ics .* cos.(op[node, :φ])
    # PD_ics[2, :] = raw_ics .* sin.(op[node, :φ])
    # base_state = copy(op.vec)
    # base_state[1:2] .= 0.
    # ics = ICset(base_state, PD_ics, 1:2)

    # @time (μ, μerr, total)= basin_stability_fixpoint(
    #     pg,
    #     op.vec,
    #     ics;
    #     distance = pg_dist,
    #     threshold = 1E-4,
    #     parallel_alg = EnsembleThreads(),
    #     solver = Rodas5(),
    #     verbose = true,
    #     return_df = false,
    #     )

    # amp = μ ± μerr

    # # perturb phase

    # lb = op[node, :φ] - perturbations[:φ]
    # ub = op[node, :φ] + perturbations[:φ]

    # raw_ics = QuasiMonteCarlo.sample(sample_size, lb, ub, SobolSample())

    # PD_ics = zeros(2, sample_size)
    # PD_ics[1, :] = op[node, :v] .* cos.(raw_ics)
    # PD_ics[2, :] = op[node, :v] .* sin.(raw_ics)
    # base_state = copy(op.vec)
    # base_state[1:2] .= 0.
    # ics = ICset(base_state, PD_ics, 1:2)

    # @time (μ, μerr, total)= basin_stability_fixpoint(
    #     pg,
    #     op.vec,
    #     ics;
    #     distance = pg_dist,
    #     threshold = 1E-4,
    #     parallel_alg = EnsembleThreads(),
    #     solver = Rodas5(),
    #     verbose = true,
    #     return_df = false,
    #     )

    # phase = μ ± μerr

    return combined # , freq, amp, phase
end

# very long expression ends here
########################################################
##

sample_size = 1000

node = "bus1" # schmietendorf
model = "mixed_IEEE_variation" * "_" * node

vars = [:φ, :v, :ω]
dimensions = 1:3

perturbations = Dict(:φ => π, :ω => 10, :v => 1)

sf_interval = range(0, 2, length=11)

#amp = zeros(Measurement, length(sf_interval))
#phase = zeros(Measurement, length(sf_interval))
#freq = zeros(Measurement, length(sf_interval))
combined = zeros(Measurement, length(sf_interval))

#amp_approx = zeros(Measurement, length(sf_interval))
#phase_approx = zeros(Measurement, length(sf_interval))
#freq_approx = zeros(Measurement, length(sf_interval))
combined_approx = zeros(Measurement, length(sf_interval))
##

# the core loop

for (run, var) in enumerate(sf_interval)
    @show (run, var)

    @show stability_interval(var)

    pg, op, pg_dist = original_model(;scaling=var)

    combined[run] = dimensional_basin_analyses(pg, op, pg_dist)

    pg, op, pg_dist = approximate_model(;scaling=var)

    combined_approx[run] = dimensional_basin_analyses(pg, op, pg_dist)
    #, freq_approx[run], amp_approx[run], phase_approx[run]
end


#df = DataFrame(c = combined, f=freq, a=amp, p=phase, ca = combined_approx, fa=freq_approx, aa=amp_approx, pa=phase_approx)
df = DataFrame(c = combined, ca = combined_approx)
CSV.write("$(plot_path)$(model)_scaling.csv", df)

##
df = CSV.read("$(plot_path)$(model)_scaling.csv", DataFrame)
c = parse.(Measurement{Float64},df.c)
ca = parse.(Measurement{Float64},df.ca)

min_v = map(x -> stability_interval(x)[1:2], sf_interval)

# simple plot
p = plot(sf_interval, c, c=:black, xguide=L"$f_S$", yguide=L"\mu_1", ylims=(0, 1), label="original system", legend=:bottomleft, right_margin=15Plots.mm)
plot!(p, sf_interval,ca, c=:black, ls=:dash, label="normal form")
p2 = twinx()
plot!(p2, sf_interval, first.(min_v), yguide=L"\rho_{min}\;[pu]", ylim=(0,1), c=:red, label=false, guidefontcolor=:red)
plot!(p2, sf_interval, last.(min_v), c=:red, ls=:dash, label=false)
plot!(tickfont=12, 
    guidefont=14, 
    legendfont=12)
# plot!(sf_interval, df.f, c=:blue, label="ω")
# plot!(sf_interval, df.fa, c=:blue, ls=:dash, label=false)
# plot!(sf_interval, df.a, c=:red, label="ρ")
# plot!(sf_interval, df.aa, c=:red, ls=:dash, label=false)
# plot!(sf_interval, df.p, c=:orange, label="φ")
# plot!(sf_interval, df.pa, c=:orange, ls=:dash, label=false)


savefig("$(plot_path)$(model)_sf.png")
