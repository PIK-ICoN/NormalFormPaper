## This is the script for the streamplots in the paper.

using GLMakie
using LaTeXStrings

##

plot_path = "$(@__DIR__)/../figures/"
fig = Figure(resolution=(1500, 700),fontsize=30, font="CMU Serif")

## normal form parameters (no internal variable)

Ar = 1.0; Ai = 1.0
Cr = -1.0; Ci = 0.0
Gr = 0.0; Gi = 1.0
Hr = 1.0; Hi = 0.0

## normal form model

model(u,δs) = Point2f(
    Ar*u[1] - Ai*u[2] + Cr*(u[1]^2 + u[2]^2)*u[1] - Ci*(u[1]^2 + u[2]^2)*u[2] + Gr*δs[1]*u[1] - Gi*δs[1]*u[2] + Hr*δs[2]*u[1] - Hi*δs[2]*u[2],
    Ar*u[2] + Ai*u[1] + Cr*(u[1]^2 + u[2]^2)*u[2] + Ci*(u[1]^2 + u[2]^2)*u[1] + Gr*δs[1]*u[2] + Gi*δs[1]*u[1] + Hr*δs[2]*u[2] + Hi*δs[2]*u[1],
)

## first plot

δp = 0.
δq = 0.
obs = Base.Fix2(model,[δp,δq])

ax1 = Axis(fig[1, 1], xlabel = L"\Re(u) \;(pu)", ylabel = L"\Im(u) \;(pu)",aspect = 1)
sp = streamplot!(obs,
                -1.5..1.5, -1.5..1.5,
                gridsize = (20, 20),
                arrow_size = 0.07,
                linewidth = 2.5,
                colormap = :inferno,
                colorrange = (0,6),
                gridalpha=0.)
hidedecorations!(ax1, ticks=false, ticklabels = false, label = false)

## second plot

δp = -0.9
δq = 0.
obs = Base.Fix2(model,[δp,δq])

ax2 = Axis(fig[1, 2], xlabel = L"\Re(u) \;(pu)", ylabel = L"\Im(u) \;(pu)", aspect = 1)
sp = streamplot!(obs,
                -1.5..1.5, -1.5..1.5,
                gridsize = (20, 20),
                arrow_size = 0.07,
                linewidth = 2.5,
                colormap = :inferno,
                colorrange = (0,6))
hidedecorations!(ax2, ticks=false, ticklabels = false, label = false)

## third plot

δp = 0
δq = 0.7
obs = Base.Fix2(model,[δp,δq])

ax3 = Axis(fig[1, 3], xlabel = L"\Re(u) \;(pu)", ylabel = L"\Im(u) \;(pu)", aspect = 1)
sp = streamplot!(obs,
                -1.5..1.5, -1.5..1.5,
                gridsize = (20, 20),
                arrow_size = 0.7,
                linewidth = 2.5,
                colormap = :inferno,
                colorrange = (0,6))
hidedecorations!(ax3, ticks=false, ticklabels = false, label = false)

## add colorbar

Colorbar(fig[2, 1:3],
        limits = (0,6),
        colormap=:inferno,
        vertical=false,
        flipaxis = false,
        width = 800,
        size = 20,
        label=L"\frac{du}{dt} \; \left(\frac{pu}{s}\right)")

save(plot_path * "streamplot.png",fig)