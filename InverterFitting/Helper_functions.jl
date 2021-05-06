using LinearAlgebra, DelimitedFiles, FFTW
using Interpolations

include("Definitions.jl")

struct measurement_data
    t
    Δt
    tspan
    tsteps
    ur
    ui
    ir_fun
    ii_fun
    ω
end

function parse_data(path,tmax;detrend=0.0,angle=0.0)

  ### read csv files and transform to phasor representation ###
  file = collect(eachrow(readdlm(path, ',', Float64, header=false)))
  fac = 1.087/sqrt(3.0)
  inverter_meas = measurement(file[1], fac*file[14:15]..., -fac*(file[14].+file[15]),
                                -(file[5:7].+file[11:13])...)
  inverter = toComplex(inverter_meas)

  ### get measurement arrays (use co-rotating frame, cut off, average & change resolution) ###
  nsteps = tmax*1e3
  t = inverter.t[1:100:Int(nsteps)]
  Δt = t[2]
  tspan = (t[1],t[end])
  tsteps = t[1]:Δt:t[end]
  
  ur = averageWin(real(inverter.uRot), 100)[1:100:Int(nsteps)]
  ui = averageWin(imag(inverter.uRot), 100)[1:100:Int(nsteps)]
  ir = averageWin(real(inverter.iRot), 100)[1:100:Int(nsteps)]
  ii = averageWin(imag(inverter.iRot), 100)[1:100:Int(nsteps)]

  ### detrend arrays ###
  ur, ui = detrend_array(ur,ui,t,detrend,angle)
  ir, ii = detrend_array(ir,ii,t,detrend,angle)

  ### calculate frequency ###
  ϕ = getAngle(ur+1im*ui)
  ω = diff(vcat(ϕ[1],ϕ))/Δt

  # ρ = sqrt.(ur.^2 + ui.^2)
  # P = ur.*ir + ui.*ii
  # Q = ui.*ir - ur.*ii

  ### input functions (interpolate measurement array) ###
  itp_ir = interpolate(ir, BSpline(Linear()))
  etp_ir = extrapolate(itp_ir, Interpolations.Flat()) #extrapolate first and last value to all times
  itp_ii = interpolate(ii, BSpline(Linear()))
  etp_ii = extrapolate(itp_ii, Interpolations.Flat())
  ir_fun(t) = etp_ir(t/Δt .+1)
  ii_fun(t) = etp_ii(t/Δt .+1)

  measurement_data(t,Δt,tspan,tsteps,ur,ui,ir_fun,ii_fun,ω)

end

function detrend_array(xr,xi,t,ω,δ)
  yr = cos.(ω*t.+δ).*xr - sin.(ω*t.+δ).*xi
  yi = sin.(ω*t.+δ).*xr + cos.(ω*t.+δ).*xi
  return yr, yi
end

function combine_data(meas1,meas2,tstart1,tstart2,tend1,tend2)
    
    @assert meas1.t == meas2.t
    Δt = meas1.Δt
    tlen1 = tend1 - tstart1
    tlen2 = tend2 - tstart2
    tlen = tlen1 + tlen2 - Δt
    t = collect(0:Δt:tlen)
    tspan = (t[1],t[end])
    tsteps = t[1]:Δt:t[end]

    s1 = Int(tstart1/Δt)+1
    s2 = Int(tstart2/Δt)+1
    e1 = Int(tend1/Δt)
    e2 = Int(tend2/Δt)

    ur = vcat(meas1.ur[s1:e1],meas2.ur[s2:e2])
    ui = vcat(meas1.ui[s1:e1],meas2.ui[s2:e2])
    ω = vcat(meas1.ω[s1:e1],meas2.ω[s2:e2])

    ir_fun(t) = t < tlen1 ? meas1.ir_fun(t+tstart1) : meas2.ir_fun(t-tlen1+tstart2)
    ii_fun(t) = t < tlen1 ? meas1.ii_fun(t+tstart1) : meas2.ii_fun(t-tlen1+tstart2)

    measurement_data(t,Δt,tspan,tsteps,ur,ui,ir_fun,ii_fun,ω)

end

function plotio(meas)
    t = meas.t
    ir = meas.ir_fun.(t)
    ii = meas.ii_fun.(t)
    ur = meas.ur
    ui = meas.ui
    ω = meas.ω
    plt1 = plot(t,ir); plot!(plt1,t,ii);
    plt2 = plot(t,ur); plot!(plt2,t,ui);
    plt3 = plot(t,ω)
    plot(plt1,plt2,plt3,layout = (3,1),size=(900,900),legend=false) |> display
end

function plotPQrf(meas)
  t = meas.t
  ir = meas.ir_fun.(t); ii = meas.ii_fun.(t)
  i = ir + 1im*ii
  ur = meas.ur; ui = meas.ui
  u = ur + 1im*ui
  ω = meas.ω
  ρ = sqrt.(abs.(u.*conj.(u)))
  S = u.*conj.(i)
  P = real.(S); Q = imag.(S)
  plt1 = plot(t,P)
  plt2 = plot(t,Q)
  plt3 = plot(t,ω)
  plt4 = plot(t,ρ)
  plot(plt1,plt2,plt3,plt4,layout = (4,1),size=(900,900),legend=false) |> display
end


function get_variables(sol::ODESolution)
    ur = sol[1,:]
    ui = sol[2,:]
    ω = sol[3,:]
    ρ = sqrt.(ur.^2+ui.^2)
    t = sol.t
    return ur, ui, ω, ρ, t
end

function droop_model_pars(τp,kp,kq,Vd,V0,P0,Q0,ω0)
    P = zeros(15)
    P[1] = V0^2/(2*τp*Vd^2) + kq*Q0/(τp*Vd)
    P[2] = -ω0
    P[3] = 0
    P[4] = 1
    P[5] = -1/(2*τp*Vd^2)
    P[6] = 0
    P[7] = 0
    P[8] = 0
    P[9] = -kq/(τp*Vd)
    P[10] = 0
    P[11] = ω0/τp + kp*P0/τp
    P[12] = -1/τp
    P[13] = 0
    P[14] = -kp/τp
    P[15] = 0
    return P
end

function parameter_mask(P)
  P[2] = 0
  P[3] = 0
  #P[4] = 1
  #P[5] = -1/(2*τ)
  P[6] = 0
  P[7] = 0
  P[8] = 0
  P[10] = 0
  #P[12] = P[12] = -1/τ
  P[13] = 0
  P[15] = 0
  return P
end
    