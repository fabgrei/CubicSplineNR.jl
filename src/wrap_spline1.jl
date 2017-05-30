@compat abstract type CubicSpline{T <: Real} end

immutable ImmutableCubicSpline{T<:Real} <: CubicSpline{T}
    x::Vector{T}
    fx::Vector{T}
    fdp::Vector{T}
    npts::Ref{Cint}
end

type MutableCubicSpline{T<:Real} <: CubicSpline{T}
    x::Vector{T}
    fx::Vector{T}
    fdp::Vector{T}
    npts::Ref{Cint}
    dirty::Bool
end

const spline_lib = Pkg.dir()*"\CubicSplineNR\lib\spline1"

"""
    _dodp!(...)

    calls Fortran function from Numerical Recipes, provided by Tony
    Smith

    It computes the second derivatives needed for calculating the
    cubic spline interpolation.
"""
function _dodp!(x, fx, fdp, npts)
    ccall((:__procedures_MOD_dodp, spline_lib),
       Void,
       (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}),
       x, fx, fdp, npts)
end


function CubicSpline(x::Vector{Float64}, fx::Vector{Float64}; mutable::Bool=false)
    N = length(x)
    npts = Ref{Cint}(N)
    fdp = zeros(N)

    _dodp!(x, fx, fdp, npts)

    if mutable
        return MutableCubicSpline(copy(x), copy(fx), fdp, npts, false)
    else
        return ImmutableCubicSpline(x, fx, fdp, npts)
    end
end

function update!(cs::MutableCubicSpline, fx::Vector{Float64})
    cs.fdp .= 0.0
    cs.fx .= fx .+ 0.0

    dodp!(cs.x, cs.fx, cs.fdp, cs.npts)

    cs.dirty = false
end

function update!(cs::MutableCubicSpline, x::Vector{Float64}, fx::Vector{Float64})
    cs.fdp .= 0.0
    cs.x .= x .+ 0.0
    cs.fx .= fx .+ 0.0

    dodp!(cs.x, cs.fx, cs.fdp, cs.npts)

    cs.dirty = false
end

# function set!(cs::MutableCubicSpline, i::Int, value::Float64)
#     cs.fx[i] = value
#     cs.dirty = true
# end

function _interp!(point::Float64,
                 x::Vector, fx::Vector, fdp::Vector, npts,
                 level::Int, deriv::Int, second::Int,
                 y::Vector{Float64}, yp::Vector{Float64}, ydp::Vector{Float64})
    ny = Ref{Cint}(level)
    nyp = Ref{Cint}(deriv)
    nydp = Ref{Cint}(second)

    point = Ref{Cdouble}(point)

    ccall((:__procedures_MOD_interp, "spline1"), Void,
(Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
    point, x, fx, fdp, y, yp, ydp, ny, nyp, nydp, npts)
end


function interp_level{CS <: CubicSpline{Float64}}(cs::CS, x::Float64)
    point = Ref{Cdouble}(x)

    y = zeros(1)
    yp = zeros(1)
    ydp = zeros(1)

    _interp!(x, cs.x, cs.fx, cs.fdp, cs.npts, 1, 0, 0, y, yp, ydp)
#     ccall((:__procedures_MOD_interp, "spline1"), Void,
# (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
#     point, cs.x, cs.fx, cs.fdp, y, yp, ydp, ny, nyp, nydp, cs.npts)
     y[1]
end

#function interp_level{T<:Float64, CS <: CubicSpline{T}}(cs::CS, x::ForwardDiff.Dual{T})
function interp_level{T<:Float64, CS <: CubicSpline{Float64}}(cs::CS, x::ForwardDiff.Dual{T})

    ny = Ref{Cint}(1)
    nyp = Ref{Cint}(1)
    nydp = Ref{Cint}(0)

    point = Ref{Cdouble}(ForwardDiff.value(x))

    y = zeros(1)
    yp = zeros(1)
    ydp = zeros(1)

    ccall((:__procedures_MOD_interp, "spline1"), Void,
(Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
    point, cs.x, cs.fx, cs.fdp, y, yp, ydp, ny, nyp, nydp, cs.npts)
    ForwardDiff.Dual(y[1]::T, yp[1]::T * ForwardDiff.partials(x))
end

function(itp::ImmutableCubicSpline)(x::Real) #where T <: Real
    interp_level(itp, x)
end

function(itp::MutableCubicSpline)(x::Real; update::Bool=true) #where T <: Real
    if update && itp.dirty
        update!(itp, itp.fx)
    end
    interp_level(itp, x)
end
