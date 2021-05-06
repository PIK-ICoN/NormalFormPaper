##

using OrdinaryDiffEq
using Statistics
distance(x,y) = sqrt(mean((x .- y) .^ 2))
using Flux, DiffEqFlux
using Plots

##

include("Helper_functions.jl")
path = ("./Data/")
tmax = 70 #seconds
meas1 = parse_data(joinpath(path,"Testcase5/TC5.13.txt"),tmax,detrend=0.062)
meas2 = parse_data(joinpath(path,"Testcase5/TC5.23.txt"),tmax,detrend=-0.00245,angle=1.816-π/2)#angle=-0.065)
meas3 = parse_data(joinpath(path,"Testcase5/TC5.11.txt"),tmax,detrend=0.051,angle=-π/2)
##

### combine measurements ###
#meas = combine_data(meas3,meas2,0,0,55,25)
meas = meas3
#plotio(meas)
#plotPQrf(meas)

##

u_l = meas.ur .+ 1.0im .* meas.ui
i_l = meas.ir_fun.(meas.t) + 1.0im * meas.ii_fun.(meas.t)
ω_l = meas.ω

dt = meas.Δt
du = diff(u_l)/dt
dω = diff(ω_l)/dt

u = u_l[1:end-1]
i = i_l[1:end-1]
ω = ω_l[1:end-1]

du_r = real.(du)
du_i = imag.(du)
u_r = real.(u)
u_i = imag.(u)

##

S = u .* conj.(i)
c_term = ones(length(S))
P = real.(S)
Q = imag.(S)
rho = abs.(u .* conj.(u))

##

A = [c_term ω rho P Q]

A_comp = [diagm(u_r) * A -diagm(u_i) * A zeros(size(A));
          diagm(u_i) * A  diagm(u_r) * A zeros(size(A));
          zeros(size(A)) zeros(size(A)) A]

# du_r = A * diag(u_r) * x_r .- A * diag(u_i) * x_i
# du_i = A * diag(u_i) * x_r .+ A * diag(u_r) * x_i 
# dω   = A * x_ω

##

du_c = [du_r; du_i; dω]

p_comp = A_comp \ du_c

##

du_t_c = A_comp * p_comp

@info distance(du_t_c, du_c)

##

u0 = [meas.ur[1];meas.ui[1];meas.ω[1]]
p_lin = zeros(15)
p_lin[1:2:10] .= p_comp[1:5]
p_lin[2:2:10] .= p_comp[6:10]
p_lin[11:15] .= p_comp[11:15]

function model(du,u,iᵣ,iᵢ,p,t)
    uᵣ, uᵢ, x = u
    P = uᵣ*iᵣ(t) + uᵢ*iᵢ(t)
    Q = uᵢ*iᵣ(t) - uᵣ*iᵢ(t)
    ρ2 = uᵣ^2 + uᵢ^2
    du[1] = duᵣ = p[1]*uᵣ - p[2]*uᵢ + p[3]*x*uᵣ - p[4]*x*uᵢ + p[5]*uᵣ*ρ2 - p[6]*uᵢ*ρ2 + p[7]*uᵣ*P - p[8]*uᵢ*P + p[9]*uᵣ*Q - p[10]*uᵢ*Q
    du[2] = duᵢ = p[1]*uᵢ + p[2]*uᵣ + p[3]*x*uᵢ + p[4]*x*uᵣ + p[5]*uᵢ*ρ2 + p[6]*uᵣ*ρ2 + p[7]*uᵢ*P + p[8]*uᵣ*P + p[9]*uᵢ*Q + p[10]*uᵣ*Q
    du[3] = dx = p[11] + p[12]*x + p[13]*ρ2 + p[14]*P + p[15]*Q
  end

##

prob = ODEProblem((du,u,p,t) -> model(du,u,meas.ir_fun,meas.ii_fun,p,t),u0,meas.tspan,p_lin)
sol = solve(prob, Tsit5(),saveat=meas.tsteps)

##

data = hcat(meas.ur,meas.ui,meas.ω) |> transpose
t = meas.tsteps
ρm = sqrt.(data[1,:].^2 + data[2,:].^2)
ρs = sqrt.(sol[1,:].^2 + sol[2,:].^2)
plt1 = plot(t,data[1,:]); plot!(plt1,sol,vars=1)
plt2 = plot(t,data[2,:]); plot!(plt2,sol,vars=2)
plt3 = plot(t,data[3,:]); plot!(plt3,sol,vars=3)
plt4 = plot(t,ρm); plot!(plt4,t,ρs)
plot(plt1,plt2,plt3,plt4,layout = (4,1),size=(900,900),legend=false) |> display

##

function predict(u0,p,meas)
    prob = ODEProblem((du,u,p,t) -> model(du,u,meas.ir_fun,meas.ii_fun,p,t),u0,meas.tspan,p)
    sol = solve(prob,Tsit5(),saveat=meas.tsteps)
    return sol
end
  
function loss(par)
    sol = predict(par[1:3],par[4:18],meas)
    ur, ui, ω1, ρ1, t = get_variables(sol)
    ω2 = data[3,:]
    ρ2 = sqrt.(data[1,:].^2 + data[2,:].^2)
    loss = sum(abs.(meas.ur-ur))
    loss += sum(abs.(meas.ui-ui))
    loss += sum(abs.(ω1-ω2))
    loss += 10*sum(abs.(ρ1-ρ2))
    return loss, Array(sol)
end
  
cb = function (par, loss, sol)
    display(loss)
    ρm = sqrt.(meas.ur.^2 + meas.ui.^2)
    ρs = sqrt.(sol[1,:].^2 + sol[2,:].^2)
    plt1 = plot(meas.t,meas.ur); plot!(plt1,meas.t,sol[1,:])
    plt2 = plot(meas.t,meas.ui); plot!(plt2,meas.t,sol[2,:])
    plt3 = plot(meas.t,meas.ω); plot!(plt3,meas.t,sol[3,:])
    plt4 = plot(meas.t,ρm); plot!(plt4,meas.t,ρs)
    display(plot(plt1,plt2,plt3,plt4,layout = (4,1),size=(900,900),legend=false))
    return false
end

##
  
par = zeros(18)
par[1:3] = u0; par[4:18] = p_lin;

##
res = DiffEqFlux.sciml_train(loss, par,
                                    ADAM(0.01),
                                    cb = cb,
                                    maxiters = 100)
                                
##


par = res.minimizer

##

u0 = par[1:3]; u0 = [meas.ur[1] meas.ui[1] 0]
prob = ODEProblem((du,u,p,t) -> model(du,u,meas.ir_fun,meas.ii_fun,p,t),u0,meas.tspan,par[4:18])
sol = solve(prob, Tsit5(),saveat=meas.tsteps)

##

data = hcat(meas.ur,meas.ui,meas.ω) |> transpose
t = meas.tsteps
ρm = sqrt.(data[1,:].^2 + data[2,:].^2)
ρs = sqrt.(sol[1,:].^2 + sol[2,:].^2)
plt1 = plot(t,data[1,:]); plot!(plt1,sol,vars=1); ylabel!(plt1,"Re(u)")
plt2 = plot(t,data[2,:]); plot!(plt2,sol,vars=2); ylabel!(plt2,"Im(u)")
plt3 = plot(t,data[3,:]); plot!(plt3,sol,vars=3); ylabel!(plt3,"ω")
plt4 = plot(t,ρm); plot!(plt4,t,ρs); xlabel!(plt4,"t [s]"); ylabel!(plt4,"|u|")
plot(plt1,plt2,plt3,plt4,layout = (4,1),size=(700,700),legend=false,linewidth=2) |> display

##

#savefig("modelfit.png")
#@save "min_parameters.jld2" par