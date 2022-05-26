using PowerDynamics

include("components.jl")

##
# line parameters from hannes' bachelor thesis which derives them from the Schmietendorf paper
B = 1.
Bshunt = B - 0.8
ic_guess = [1, 0, 0, 1, 0, 0]

lines = [StaticLine(from=1, to=2, Y = - B * 1im), ]

Sn = 0.505

##

#pq = PowerGrid([PQAlgebraic(;P=0.5, Q=0.2), PQAlgebraic(; P=-0.5, Q=0.2)], lines) |> find_operationpoint

##

nodes_orig = [
    SchifferOriginal(;τ_P=2.5, τ_Q=2.5, K_P=5, K_Q=0.1, V_r=1, P=0.5, Q=0.2, Y_n=Bshunt*1im), # Hannes VVogel BA
    #SchifferOriginal(;τ_P=2.5, τ_Q=2.5, K_P=5, K_Q=0.1, V_r=1, P=-0.5, Q=0.2, Y_n=Bshunt*1im),
    # SchifferOriginal(;τ_P=0.5, τ_Q=0.5, K_P=0.2/Sn, K_Q=0.1/Sn, V_r=1, P=0.3*Sn, Q=0.025*Sn), # Schiffer paper
    # SchifferOriginal(;τ_P=0.5, τ_Q=0.5, K_P=0.2/Sn, K_Q=0.1/Sn, V_r=1, P=-0.3*Sn, Q=0.025*Sn),
    SlackAlgebraic(; U=complex(1)),
]


pg_orig = PowerGrid(nodes_orig, lines)
op_orig = find_operationpoint(pg_orig; sol_method=:rootfind)

@show op_orig[:, :φ];
@show op_orig[:, :v];
@show op_orig[1, :ω];
@show op_orig[:, :p];
@show op_orig[:, :q];
@show op_orig[:, :i];

##

nodes_approx = [
    SchifferApprox(;τ_P=2.5, τ_Q=2.5, K_P=5, K_Q=0.1, V_r=1, P=0.5, Q=0.2, Y_n=Bshunt*1im), # Hannes VVogel BA
    #SchifferApprox(;τ_P=2.5, τ_Q=2.5, K_P=5, K_Q=0.1, V_r=1, P=-0.5, Q=0.2, Y_n=Bshunt*1im),
    # SchifferApprox(;τ_P=0.5, τ_Q=0.5, K_P=0.2/Sn, K_Q=0.1/Sn, V_r=1, P=0.3*Sn, Q=0.025*Sn), # Schiffer paper
    # SchifferApprox(;τ_P=0.5, τ_Q=0.5, K_P=0.2/Sn, K_Q=0.1/Sn, V_r=1, P=-0.3*Sn, Q=0.025*Sn),
    SlackAlgebraic(; U=complex(1)),
]

pg_approx = PowerGrid(nodes_approx, lines)
op_approx = find_operationpoint(pg_approx; sol_method=:rootfind)

@show op_approx[:, :φ];
@show op_approx[:, :v];
@show op_approx[:, :p];
@show op_approx[:, :q];
@show op_approx[:, :i];

##
# The third-order model is basically already in the normal form, hence we don't see any differences ...


using Plots
using LaTeXStrings

plot_path = "$(@__DIR__)/../figures/"

##

tspan = (0., 50.)
tarray = 0:0.01:last(tspan)
pp = PowerPerturbation(node=1, fault_power=1, tspan_fault=(10., 12.), var=:(P))

sol_orig = simulate(pp, pg_orig, op_orig, tspan);
plot(sol_orig.(tarray, 1, :u_r), sol_orig.(tarray, 1, :u_i), label="original model")

sol_approx = simulate(pp, pg_approx, op_approx, tspan);
plot!(sol_approx.(tarray, 1, :u_r), sol_approx.(tarray, 1, :u_i), label="approximation")

scatter!(sol_orig.([pp.tspan_fault...], 1, :u_r), sol_orig.([pp.tspan_fault...], 1, :u_i), label=false, c=:black)
scatter!(sol_approx.([pp.tspan_fault...], 1, :u_r), sol_approx.([pp.tspan_fault...], 1, :u_i), label=false, c=:black)

plot!(op_orig[1, :v] .* cos.(0:0.1:π/2), op_orig[1, :v] .* sin.(0:0.1:π/2), c=:black, ls=:dash, label=false, xguide=L"$\Re (u)$", yguide=L"$\Im (u)$", legend=:bottom)

savefig(plot_path * "schiffer_traj_ru_iu.png")

##
symmod(x) = mod2pi(x + π) -π

pv = plot(tarray, t -> sol_orig.(t, 1, :v), legend=:bottomright, label="original model", ylabel=L"\rho \;(pu)", yshowaxis=true, xticks=false, yticks=[0.98, 1., 1.02]);
plot!(pv, tarray, t -> sol_approx.(t, 1, :v), label="normal form", linestyle=:dash);
vline!(pv, [pp.tspan_fault...], c=:black, alpha=0.8, label="");

pφ = plot(tarray, t -> sol_orig.(t, 1, :φ) .|> symmod, legend=false, ylabel=L"\phi \;(rad)", yshowaxis=true, xticks=false);
plot!(pφ, tarray, t -> sol_approx.(t, 1, :φ) .|> symmod, linestyle=:dash);
vline!(pφ, [pp.tspan_fault...], c=:black, alpha=0.8);

pω = plot(tarray, t -> sol_orig.(t, 1, :ω), legend=false, xlabel=L"t \;(s)", ylabel=L"\omega\; (rad/s)", yticks=[-0.5, 0., 0.5]);
plot!(pω, tarray, t -> sol_approx.(t, 1, :ω), linestyle=:dash);
vline!(pω, [pp.tspan_fault...], c=:black, alpha=0.8);

l = @layout [a{0.32h}; b{0.32h}; c{0.32h}]

plot(pv, pφ, pω;
    layout=l,
    size=(600, 600),
    grid=:y, gridalpha=0.,
    link=:x,
    fontfamily="Computer Modern",
    tickfontsize=24,
    guidefontsize=24,
    legendfontsize=18,
    linewidth=4,
)

savefig(plot_path * "schiffer_traj.png")

##

pv = plot(tarray, t -> 100 .- 100 .* sol_approx.(t, 1, :v) ./ sol_orig.(t, 1, :v), ylabel="\\Delta\\rho [%]", yshowaxis=true, xticks=false);
vline!(pv, [pp.tspan_fault...], c=:black, alpha=0.8, label="");

pφ = plot(tarray, t -> 100 .- 100 .* sol_approx.(t, 1, :φ) ./ sol_orig.(t, 1, :φ), ylabel="\\Delta\\varphi [%]", yshowaxis=true, xticks=false);
vline!(pφ, [pp.tspan_fault...], c=:black, alpha=0.8, label="");

pω = plot(tarray, t -> 100 .- 100 .* (50*2π .+ sol_approx.(t, 1, :ω)) ./ (50*2π .+ sol_orig.(t, 1, :ω)), xlabel="t", ylabel="\\Delta\\omega [%]");
vline!(pω, [pp.tspan_fault...], c=:black, alpha=0.8, label="");

l = @layout [a; b; c]
plot(pv, pφ, pω; layout=l, size=(800, 400), legend=false, grid=:y, gridalpha=0.5, link=:x)

savefig(plot_path * "schiffer_traj_error.png")


##

# git clone git@github.com:luap-pik/ProbabilisticStability.git
# ]dev path/to/repo

using Distances

# We identify convergence given the frequency vanishes and ignore phase shifts.

# since the second node is a slack we only compare the first one

function pg_dist_orig(x, y; d = Euclidean())
    sx = State(pg_orig, x)
    sy = State(pg_orig, y)
    return d([sx[1, :v]; sx[1, :ω]], [sy[1, :v]; sy[1, :ω]])
end

function pg_dist_approx(x, y; d = Euclidean())
    sx = State(pg_approx, x)
    sy = State(pg_approx, y)
    return d([sx[1, :v]; sx[1, :ω]], [sy[1, :v]; sy[1, :ω]])
end

sample_size = 40000 # square number!!

model = "schiffer"

##

# attention: delete csv files when you changed parameters in the model
include("basin_runs.jl")
