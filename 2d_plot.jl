using WGLMakie
WGLMakie.activate!()
##

f(x) = Point2f0(
    x[2] - x[1] * (1 - (x[1]^2 + x[2]^2)),
    -x[1] - x[2] * (1 - (x[1]^2 + x[2]^2)),
)

##

streamplot(f, -2.5..2.5, -2.5..2.5, colormap = :magma, gridsize = (100, 100), arrow_size = 0.01)

