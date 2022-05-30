using Interpolations
using LaTeXStrings

function time_interpolation(i,Δt)
  itp_i = interpolate(i, BSpline(Linear()))
  etp_i = extrapolate(itp_i, Interpolations.Flat()) #extrapolate first and last value to all times
  t -> etp_i(t/Δt .+1)
end

function plot_io(meas)
  t = meas.t
  ir = meas.ir
  ii = meas.ii
  ur = meas.ur
  ui = meas.ui
  ω = meas.omega
  plt1 = plot(t,ir,label = L"\Re(j)"); plot!(plt1,t,ii,label = L"\Im(j)"); ylabel!(plt1,L"j \;(pu)"), xlabel!(plt1,""); yticks!(plt1,-1.5:1.5:3.)
  plt2 = plot(t,ur,label = L"\Re(u)"); plot!(plt2,t,ui,label = L"\Im(u)"); ylabel!(plt2,L"u \;(pu)"), xlabel!(plt2,L"t \;(s)")
  plot(plt1,plt2,
      layout = (2,1),
      size=(600,300),
      legend=true,
      linewidth=2,
      fontfamily="Computer Modern",
      gridalpha=0.,
      tickfont=14, 
      guidefont=16, 
      legendfont=14,
      bottom_margin=4Plots.mm) |> display
end

function plot_sim(meas,sol)
    data = [meas.ur meas.ui meas.omega] |> transpose
    t = meas.t
    ρm = sqrt.(data[1,:].^2 + data[2,:].^2)
    ρs = sqrt.(sol[1,:].^2 + sol[2,:].^2)
    plt1 = plot(t,data[1,:]); plot!(plt1,sol,vars=1); ylabel!(plt1,L"\Re(u) \;(pu)"); xlabel!(plt1,"")
    plt2 = plot(t,data[2,:]); plot!(plt2,sol,vars=2); ylabel!(plt2,L"\Im(u) \;(pu)"); xlabel!(plt2,"")
    plt3 = plot(t,data[3,:]); plot!(plt3,sol,vars=3); ylabel!(plt3,L"\omega \;(rad/s)"); xlabel!(plt3,"")
    plt4 = plot(t,ρm); plot!(plt4,t,ρs); xlabel!(plt4,L"t \;(s)"); ylabel!(plt4,L"|u| \;(pu)"); #yticks!(plt4,0.94:0.02:1.00)
    plot(plt1,plt2,plt3,plt4,
        layout = (4,1),
        size=(600,600),
        legend=false,
        linewidth=2,
        tickfont=14,
        guidefont=16,
        fontfamily="Computer Modern",
        gridalpha=0.) |> display  
end

function predict(u0,p,meas)
    dt = meas.t[2] - meas.t[1]
    tspan = (meas.t[1],meas.t[end])
    tsteps = tspan[1]:dt:tspan[end]
    prob = ODEProblem((du,u,p,t) -> dynamic_model(du,u,ir_fun,ii_fun,p,t),u0,tspan,p)
    sol = solve(prob,Tsit5(),saveat=tsteps)
    return sol
end

function get_variables(sol::ODESolution)
    ur = sol[1,:]
    ui = sol[2,:]
    ω = sol[3,:]
    ρ = sqrt.(ur.^2+ui.^2)
    t = sol.t
    return ur, ui, ω, ρ, t
end

function loss(par)
  sol = predict(par[1:3],par[4:18],meas)
  ur, ui, ω1, ρ1, t = get_variables(sol)
  data = [meas.ur meas.ui meas.omega] |> transpose
  ω2 = data[3,:]
  ρ2 = sqrt.(data[1,:].^2 + data[2,:].^2)
  loss = sum(abs.(meas.ur-ur))
  loss += sum(abs.(meas.ui-ui))
  loss += sum(abs.(ω1-ω2))
  loss += 10*sum(abs.(ρ1-ρ2))
  return loss, sol #Array(sol)
end

cb = function (par, loss, sol)
  display(loss)
  plot_sim(meas,sol)
  return false
end