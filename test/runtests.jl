using CubicSplineNR
using Base.Test

x = collect(linspace(0, 2pi, 30))
fx = sin.(x)
fx2 = cos.(x)

sin_itp = CubicSpline(x, fx)
sin_mut = CubicSpline(x, fx, mutable=true)

sin_itp(3.14)
sin_mut(3.14)

@time update!(sin_mut, fx2)

sin_mut(3.14)
