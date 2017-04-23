module CubicSplineNR

using ForwardDiff

include("wrap_spline1.jl")

export CubicSpline, MutableCubicSpline, ImmutableCubicSpline,
    interp_level, update!

# package code goes here

end # module
