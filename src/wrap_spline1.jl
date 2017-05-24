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

const spline_lib = Pkg.dir()*"/CubicSplineNR/lib/spline1"

function CubicSpline(x::Vector{Float64}, fx::Vector{Float64}; mutable::Bool=false)
    N = length(x)
    npts = Ref{Cint}(N)
    fdp = zeros(N)
   
    ccall((:__procedures_MOD_dodp, spline_lib),
       Void,
       (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}),
       x, fx, fdp, npts)
    if mutable
        return MutableCubicSpline(copy(x), copy(fx), fdp, npts, false)
    else
        return ImmutableCubicSpline(x, fx, fdp, npts)
    end
end

function update!(cs::MutableCubicSpline, fx::Vector{Float64})
    cs.fdp .= 0.0
    cs.fx .= fx .+ 0.0
    ccall((:__procedures_MOD_dodp, spline_lib),
       Void,
       (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}),
          cs.x, cs.fx, cs.fdp, cs.npts)
    cs.dirty = false
end

function update!(cs::MutableCubicSpline, x::Vector{Float64}, fx::Vector{Float64})
    cs.fdp .= 0.0
    cs.x .= x .+ 0.0
    cs.fx .= fx .+ 0.0
    ccall((:__procedures_MOD_dodp, spline_lib),
       Void,
       (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}),
          cs.x, cs.fx, cs.fdp, cs.npts)
    cs.dirty = false
end

# function set!(cs::MutableCubicSpline, i::Int, value::Float64)
#     cs.fx[i] = value
#     cs.dirty = true
# end

function interp_level{CS <: CubicSpline{Float64}}(cs::CS, x::Float64)
    ny = Ref{Cint}(1)
    nyp = Ref{Cint}(0)
    nydp = Ref{Cint}(0)

    point = Ref{Cdouble}(x)

    y = zeros(1)
    yp = zeros(1)
    ydp = zeros(1)
    
    ccall((:__procedures_MOD_interp, "spline1"), Void,
(Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
    point, cs.x, cs.fx, cs.fdp, y, yp, ydp, ny, nyp, nydp, cs.npts)
    y[1]::Float64
end

function interp_level{T<:Float64, CS <: CubicSpline{T}}(cs::CS, x::ForwardDiff.Dual{T})

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

function(itp::ImmutableCubicSpline{T<:Real})(x::T) #where T <: Real
    interp_level(itp, x)
end

function(itp::MutableCubicSpline{T<:Real})(x::T; update::Bool=true) #where T <: Real
    if update && itp.dirty
        update!(itp, itp.fx)
    end
    interp_level(itp, x)
end
