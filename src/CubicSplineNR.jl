module CubicSplineNR

using Compat
using ForwardDiff

include("wrap_spline1.jl")
include("multi-splines.jl")

export CubicSpline, MutableCubicSpline, ImmutableCubicSpline,
    interp_level, update!
export MultiCubicSpline

# package code goes here

end # module
