##
using NetworkDynamics: ODEVertex
import PowerDynamics: dimension, symbolsof, construct_vertex 
import PowerDynamics: showdefinition

# Eq 17 in the paper
@DynamicNode dVOC(Pset, Qset, Vset, Ωset, η, α, κ) begin
    @assert η > 0 "positive control parameter"
    @assert α > 0 "positive control parameter"
    @assert isreal(Ωset) "nominal frequency, i.e. 2π50 Hz"
    @assert isreal(Pset)  "real power set point"
    @assert isreal(Qset) "reactive power set point"
    @assert Vset > 0 "voltage set point"
    @assert isreal(κ) "unifrom complex phase"
    Y =  complex(Pset, -Qset) / Vset^2 # K = exp(im κ) Y
end [] begin
    #s = u * conj(i)
    du = complex(η * α, Ωset) * u - α * η * u * abs2(u) / Vset^2 + η * exp(κ * 1im) * (Y * u - i)
end

##

@DynamicNode dVOCapprox(Pset, Qset, Vset, Ωset, η, α, κ) begin
    @assert η > 0 "positive control parameter"
    @assert α > 0 "positive control parameter"
    @assert isreal(Ωset) "nominal frequency, i.e. 2π50 Hz"
    @assert isreal(Pset)  "real power set point"
    @assert isreal(Qset) "reactive power set point"
    @assert Vset > 0 "voltage set point"
    @assert isreal(κ) "unifrom complex phase"
    Cᵤ = - α * η / (Vset^2) + η * exp(κ * 1im) * complex(Pset, -Qset) / (Vset^4)
    Gᵤ = - η * exp(κ * 1im) / (Vset^2)
    Hᵤ = 1im * η * exp(κ * 1im) / (Vset^2)
    Aᵤ = 1im * Ωset  - Cᵤ * (Vset)^2 - Gᵤ * Pset - Hᵤ * Qset
end [] begin
    s = u * conj(i)
    v2 = abs2(u)
    du = ( Aᵤ + Cᵤ * v2 + Gᵤ * real(s) + Hᵤ * imag(s) ) * u
end

##

const REGULARISATION = 1000

@DynamicNode ImprovedSwingOriginal(P, A, M, Ω) begin
end [[ω, dω]] begin
    p = u * conj(i) |> real
    v = abs(u)
    dω = ( - A * ω + (P - p) * Ω / (Ω + ω) ) / M
    dv = REGULARISATION * (1. - v)
    du = (dv / v + 1im * ω ) * u
end

##

@DynamicNode ImprovedSwingApprox(P, A, M, Ω) begin
    Aᵤ = REGULARISATION
    Bᵤ = 1im 
    Cᵤ = -REGULARISATION
    Gᵤ = 0
    Hᵤ = 0
    Mₓ = M
    Aₓ = P 
    Bₓ = - A 
    Cₓ = 0
    Gₓ = - 1 
    Hₓ = 0
    @assert isreal(Mₓ) && Mₓ > 0
    @assert isreal(Aₓ)
    @assert isreal(Bₓ)
    @assert isreal(Cₓ)
    @assert isreal(Gₓ)
    @assert isreal(Hₓ)
end [[ω, dω]] begin
    s = u * conj(i)
    v2 = abs2(u)
    dω = ( Aₓ + Bₓ * ω + Cₓ * v2 + Gₓ * real(s) + Hₓ * imag(s) ) / Mₓ
    du = ( Aᵤ + Bᵤ * ω + Cᵤ * v2 + Gᵤ * real(s) + Hᵤ * imag(s) ) * u
end

## 

@DynamicNode iSL_dvoc(Pset, Qset, Vset, d, η, α, ω0, κ) begin
    @assert η > 0 "positive control parameter"
    @assert α > 0 "positive control parameter"
    @assert isreal(ω0) "nominal frequency, i.e. 2π50 Hz"
    @assert isreal(Pset)  "real power set point"
    @assert isreal(Qset) "reactive power set point"
    @assert Vset > 0 "voltage set point"
    @assert isreal(κ) "unifrom complex phase"
    Y =  complex(Pset, -Qset) / Vset^2 # K = exp(im κ) Y
    expkappa = exp(κ * 1im)
end [[ω, dω]] begin
    s = u * conj(i)
    v2 = abs2(u)
    dω = - d * (ω - ω0) + η * d * imag(expkappa * Y) - η * d * imag(expkappa * conj(s)) / v2
    du = complex(η * α, ω) * u + η * real(expkappa * Y) * u - α * η * u * v2 / Vset^2 - η * u * real(expkappa * conj(s)) / v2
end

##
@doc doc"""
```Julia
    PhaseAmplitudeOscillator(; Aᵤ, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Mₓ, Aₓ, Bₓ, Cₓ, Gₓ, Hₓ)
```
# Mathematical Representation
This is an implementation of Eqn. 12 in the paper, i.e. the expansion of a phase amplitude oscillator at the limit cycle:
```math
\dot{u}_n \approx \left( A_n^u + B_n^u x_n \right) u_n + C_n^u u_n \lvert u_n \rvert^2 +  \left( \Gamma_n^u p_n + \Lambda_n^u q_n \right) u_n \;,\\
\dot{x}_n \approx A_n^x + B_n^x x_n + C_n^x \lvert u_n \rvert^2 + \Gamma_n^x p_n + \Lambda_n^x q_n \;,
```
The parameters are determined as:
```math
\begin{aligned}
B_n^{u,x} &=  \frac{\partial g_n^{u,x}}{\partial x_n}(y_n^\circ) \;, \quad C_n^{u,x} &= \frac{\partial g_n^{u,x}}{\partial \lvert u_n \rvert^2}(y_n^\circ) \;,\\ 
\Gamma_n^{u,x} &=  \frac{\partial g_n^{u,x}}{\partial p_n}(y_n^\circ) \;, \quad \Lambda_n^{u,x} &= \frac{\partial g_n^{u,x}}{\partial q_n}(y_n^\circ)\;.
\end{aligned}
```
"""

##

@DynamicNode PhaseAmplitudeOscillator(Aᵤ, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Mₓ, Aₓ, Bₓ, Cₓ, Gₓ, Hₓ) begin
    @assert isreal(Mₓ) && Mₓ > 0
    @assert isreal(Aₓ)
    @assert isreal(Bₓ)
    @assert isreal(Cₓ)
    @assert isreal(Gₓ)
    @assert isreal(Hₓ)
end [[ω, dω]] begin
    s = u * conj(i)
    v2 = abs2(u)
    dω = ( Aₓ + Bₓ * ω + Cₓ * v2 + Gₓ * real(s) + Hₓ * imag(s) ) / Mₓ
    du = ( Aᵤ + Bᵤ * ω + Cᵤ * v2 + Gᵤ * real(s) + Hᵤ * imag(s) ) * u
end

##
# Section 4 in the paper
function iSL(; m, μu, μx, d,  λ, k, Pset, Qset, Vset, ωset, κ = π / 2)
    coup = exp(1im * κ) * conj(complex(Pset, Qset))
    return PhaseAmplitudeOscillator(;
        Aᵤ = μu * Vset^2 + λ * (cos(κ) * Pset + sin(κ) * Qset), #μ + ζ * λ * real(coup), 
        Bᵤ = 1im, 
        Cᵤ = - μu, # - (μ + (ζ - 1) * λ * real(coup)) / Vset^2, 
        Gᵤ = - λ * cos(κ), 
        Hᵤ = - λ * sin(κ), 
        Mₓ = m, 
        Aₓ = d * ωset + μx * Vset^2 + k * (sin(κ) * Pset - cos(κ) * Qset), #d * ωset + η * k * imag(coup), 
        Bₓ = - d, 
        Cₓ = - μx, #- (η - 1) * k * imag(coup) / Vset^2, 
        Gₓ = - k * sin(κ), 
        Hₓ = k * cos(κ)
        )
end

##

@DynamicNode iSLpar(m, μu, μx, d, λ, k, Pset, Qset, Vset, ωset, κ) begin
    coup = exp(1im * κ) * conj(complex(Pset, Qset))
end [[ω, dω]] begin
    s = u * conj(i)
    v2 = abs2(u)
    Aᵤ = μu * Vset^2 + p.λ * (cos(κ) * Pset + sin(κ) * Qset) #μ + ζ * λ * real(coup), 
    Bᵤ = 1im 
    Cᵤ = - μu # - (μ + (ζ - 1) * λ * real(coup)) / Vset^2, 
    Gᵤ = - p.λ * cos(κ) 
    Hᵤ = - p.λ * sin(κ) 
    Mₓ = p.m #m ### Parameter to be varied
    Aₓ = p.d * ωset + μx * Vset^2 + p.k * (sin(κ) * Pset - cos(κ) * Qset) #d * ωset + η * k * imag(coup), 
    Bₓ = - p.d 
    Cₓ = - μx #- (η - 1) * k * imag(coup) / Vset^2, 
    Gₓ = - p.k * sin(κ) 
    Hₓ = p.k * cos(κ)
    dω = ( Aₓ + Bₓ * ω + Cₓ * v2 + Gₓ * real(s) + Hₓ * imag(s) ) / Mₓ
    du = ( Aᵤ + Bᵤ * ω + Cᵤ * v2 + Gᵤ * real(s) + Hᵤ * imag(s) ) * u
end

##

SchifferOriginal(; kwargs...) = VSIMinimal(; kwargs...)

@DynamicNode SchifferApprox(τ_P, τ_Q, K_P, K_Q, V_r, P, Q) begin
    Aᵤ = (V_r + 2 * K_Q * Q) / (2 * τ_Q * V_r)
    Bᵤ = 1im 
    Cᵤ = - 1 / (2 * τ_Q * V_r^2)
    Gᵤ = 0
    Hᵤ = - K_Q / (τ_Q * V_r)
    Mₓ = τ_P
    Aₓ = K_P * P 
    Bₓ = - 1 
    Cₓ = 0
    Gₓ = - K_P 
    Hₓ = 0
    @assert isreal(Mₓ) && Mₓ > 0
    @assert isreal(Aₓ)
    @assert isreal(Bₓ)
    @assert isreal(Cₓ)
    @assert isreal(Gₓ)
    @assert isreal(Hₓ)
end [[ω, dω]] begin
    s = u * conj(i)
    v2 = abs2(u)
    dω = ( Aₓ + Bₓ * ω + Cₓ * v2 + Gₓ * real(s) + Hₓ * imag(s) ) / Mₓ
    du = ( Aᵤ + Bᵤ * ω + Cᵤ * v2 + Gᵤ * real(s) + Hᵤ * imag(s) ) * u
end

##


@DynamicNode SchifferVIOriginal(τ_P, τ_Q, τ_V, K_P, K_Q, V_r, P, Q) begin
end [[ω, dω], [ν, dν]] begin
    s = u * conj(i)
    ρ = abs(u)
    dω = (-ω - K_P * (real(s) - P)) / τ_P
    dν = (- (τ_V + τ_Q) * ν + V_r - ρ -  K_Q * (imag(s) - Q))
    du = (ν + 1im * ω ) * u
end

@DynamicNode SchifferVIApprox(τ_P, τ_Q, τ_V, K_P, K_Q, V_r, P, Q) begin
    Aᵤ = V_r^2 + 1 / (2 * τ_V * τ_Q) + K_Q *Q / (τ_V * τ_Q * V_r)
    Bᵤ = - (1 / τ_V) - (1 / τ_Q)  - 2V_r
    Cᵤ = - 1 / (2 * τ_V * τ_Q)
    Gᵤ = 0
    Hᵤ = - K_Q / (V_r * τ_V * τ_Q)
    Mₓ = τ_P
    Aₓ = K_P * P 
    Bₓ = - 1 
    Cₓ = 0
    Gₓ = - K_P 
    Hₓ = 0
    @assert isreal(Mₓ) && Mₓ > 0
    @assert isreal(Aₓ)
    @assert isreal(Bₓ)
    @assert isreal(Cₓ)
    @assert isreal(Gₓ)
    @assert isreal(Hₓ)
end [[ω, dω], [ν, dν]] begin
    s = u * conj(i)
    v2 = abs2(u)
    dω = ( Aₓ + Bₓ * ω + Cₓ * v2 + Gₓ * real(s) + Hₓ * imag(s) ) / Mₓ
    dν =   Aᵤ + Bᵤ * ν + Cᵤ * v2 + Gᵤ * real(s) + Hᵤ * imag(s) 
    du = (ν + 1im * ω ) * u
end

##


@DynamicNode SchmietendorfOriginal(P_m, γ, α, E_f, X) begin
    @assert γ > 0
end [[ω, dω], ] begin
    s = u * conj(i)
    e_q = abs(u)
    de_q = (E_f - e_q - X * imag(s) / e_q) / α
    dω = P_m - γ * ω - real(s)
    du = (de_q / e_q + 1im * ω ) * u
end

##

@DynamicNode SchmietendorfApprox(P_m, γ, α, E_f, X, E_c, Q_c) begin
    Aᵤ = ((3E_f / 2E_c) - 1 + X * Q_c / (E_c)^2 ) / α
    Bᵤ = 1im 
    Cᵤ = (- E_f / (2E_c^3) + X * Q_c / (E_c)^4 ) / α
    Gᵤ = 0
    Hᵤ = - X / (α * (E_c)^2)
    Mₓ = 1
    Aₓ = P_m 
    Bₓ = - γ
    Cₓ = 0. 
    Gₓ = -1
    Hₓ = 0
    @assert isreal(Mₓ) && Mₓ > 0
    @assert isreal(Aₓ)
    @assert isreal(Bₓ)
    @assert isreal(Cₓ)
    @assert isreal(Gₓ)
    @assert isreal(Hₓ)
end [[ω, dω]] begin
    s = u * conj(i)
    v2 = abs2(u)
    dω = ( Aₓ + Bₓ * ω + Cₓ * v2 + Gₓ * real(s) + Hₓ * imag(s) ) / Mₓ
    du = ( Aᵤ + Bᵤ * ω + Cᵤ * v2 + Gᵤ * real(s) + Hᵤ * imag(s) ) * u
end

##