##

using CSV, DataFrames
using OrdinaryDiffEq
using Plots
using Flux, DiffEqFlux

include("function_definitions.jl")

##

#=
Parse the data of the fist measurement from the csv-file and transform it into a DataFrame.
=#

meas = CSV.File("./Measurement1.csv") |> DataFrame
plot_io(meas)

##

#=
Calculate complex voltage/current, active/reactive power and numerical derivatives.
=#

u = meas.ur .+ 1im * meas.ui
i = meas.ir .+ 1im * meas.ii
ω = meas.omega

dt = meas.t[2] - meas.t[1]
du = diff(u)/dt
dω = diff(ω)/dt

u = u[1:end-1]
i = i[1:end-1]
ω = ω[1:end-1]

S = u .* conj.(i)
P = real.(S)
Q = imag.(S)
ρ = abs.(u .* conj.(u))

##

#=
Estimate the parameters of the normal form model with a linear regression.
For this we calculate the lhs of the ODE du/u and define a matrix with the time series' as coloumns:
A = [1 ω ρ P Q]
Then we separate real and imaginary parts and get the following equations for the parameters:
real(du/u) = A * p_r
imag(du/u) = A * p_i
dω = A * p_ω
We solve these equations by using the backslash operator "\" yielding a least squares approximate solution.
=#

c = ones(length(P))
A = [c ω ρ P Q]
p_r = A \ real.(du ./ u)
p_i = A \ imag.(du ./ u)
p_ω = A \ dω

p_reg = [p_r; p_i; p_ω]

##

#=
Define the dynamic model. We consider the currents to be input functions driving the voltage dynamics.
=#

function dynamic_model(du,u,ir,ii,p,t)
    ur, ui, x = u
    P = ur*ir(t) + ui*ii(t)
    Q = ui*ir(t) - ur*ii(t)
    ρ2 = ur^2 + ui^2
    du[1] = dur = p[1]*ur - p[2]*ui + p[3]*x*ur - p[4]*x*ui + p[5]*ur*ρ2 - p[6]*ui*ρ2 + p[7]*ur*P - p[8]*ui*P + p[9]*ur*Q - p[10]*ui*Q
    du[2] = dui = p[1]*ui + p[2]*ur + p[3]*x*ui + p[4]*x*ur + p[5]*ui*ρ2 + p[6]*ur*ρ2 + p[7]*ui*P + p[8]*ur*P + p[9]*ui*Q + p[10]*ur*Q
    du[3] = dx = p[11] + p[12]*x + p[13]*ρ2 + p[14]*P + p[15]*Q
end

##

#=
Calculate model parameters, initial conditions and interpolate the current signal.
=#

u0 = [meas.ur[1];meas.ui[1];meas.omega[1]]

p = zeros(15)
p[1:2:10] .= p_reg[1:5]
p[2:2:10] .= p_reg[6:10]
p[11:15] .= p_reg[11:15]

tspan = (meas.t[1],meas.t[end])
tsteps = tspan[1]:dt:tspan[end]

ir_fun = time_interpolation(meas.ir,dt)
ii_fun = time_interpolation(meas.ii,dt)

##

#=
Simulate the dynamics model using the parameter estimates from the linear regression.
=#

prob = ODEProblem((du,u,p,t) -> dynamic_model(du,u,ir_fun,ii_fun,p,t),u0,tspan,p)
sol = solve(prob, Tsit5(),saveat=tsteps)
plot_sim(meas,sol)

##

#=
Optimize parameters and initial conditions with a stochastic gradient descent
using the DiffEqFlux package and the ADAM solver.
=#
  
par = zeros(18)
par[1:3] = u0; par[4:18] = p;

res = DiffEqFlux.sciml_train(loss, par,
                                    ADAM(0.01),
                                    cb = cb,
                                    maxiters = 100)

par = res.minimizer
savefig("modelfit1.png")
                                
##

#=
Validate the parameters by plotting the simulations for a different data set.
=#

meas = CSV.File("./Measurement2.csv") |> DataFrame
plot_io(meas)

ir_fun = time_interpolation(meas.ir,dt)
ii_fun = time_interpolation(meas.ii,dt)

u0 = par[1:3]; u0 = [meas.ur[1] meas.ui[1] 0]
prob = ODEProblem((du,u,p,t) -> dynamic_model(du,u,ir_fun,ii_fun,p,t),u0,tspan,par[4:18])
sol = solve(prob, Tsit5(),saveat=tsteps)
plot_sim(meas,sol)
savefig("modelfit2.png")

##