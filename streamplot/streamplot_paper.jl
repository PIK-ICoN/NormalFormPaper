## This is the script for the streamplots in the paper.

using GLMakie

##

plot_path = "$(@__DIR__)/../figures/"
fig = Figure(resolution=(750, 750),fontsize=30, font="CMU Serif")
ax = Axis(fig[1, 1], xlabel = L"\Re(u)", ylabel = L"\Im(u)")

## deviations of active and reactive power

δp = 0
δq = 0.7

## normal form parameters (no internal variable)

Ar = 1.0; Ai = 1.0
Cr = -1.0; Ci = 0.0
Gr = 0.0; Gi = 1.0
Hr = 1.0; Hi = 0.0

## normal form model

model(u,δs) = Point2f0(
    Ar*u[1] - Ai*u[2] + Cr*(u[1]^2 + u[2]^2)*u[1] - Ci*(u[1]^2 + u[2]^2)*u[2] + Gr*δs[1]*u[1] - Gi*δs[1]*u[2] + Hr*δs[2]*u[1] - Hi*δs[2]*u[2],
    Ar*u[2] + Ai*u[1] + Cr*(u[1]^2 + u[2]^2)*u[2] + Ci*(u[1]^2 + u[2]^2)*u[1] + Gr*δs[1]*u[2] + Gi*δs[1]*u[1] + Hr*δs[2]*u[2] + Hi*δs[2]*u[1],
)

## fix δp and δq to the values above

obs = Base.Fix2(model,[δp,δq])

## streamplot with GLMakie

sp = streamplot!(obs,-1.5..1.5, -1.5..1.5, gridsize = (20, 20), arrow_size = 0.07, linewidth = 2.5, colormap = :magma)
hidedecorations!(ax, ticks=false, ticklabels = false, label = false)
#Colorbar(fig[2,1],vertical=false, flipaxis = false, colormap = :magma)
display(fig)
save(plot_path * "streamplot3.png",fig)