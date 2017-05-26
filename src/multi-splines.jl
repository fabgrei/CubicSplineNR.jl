type MultiCubicSpline{T<:Real}
    x::Vector{T}
    fx::Matrix{T}
    fdp::Matrix{T}
    npts::Ref{Cint}
    dirty::Bool
end

function CubicSpline(x::Vector{Float64}, fx::Matrix{Float64})
    N, M = size(fx)
    npts = Ref{Cint}(N)
    fdp = zeros(fx)

    for i in 1:M
        _dodp!(x, view(fx, :, i), view(fdp, :,i), npts)
    end

    MultiCubicSpline(x, fx, fdp, npts, false)
end

function update!(cs::MultiCubicSpline, fx::Matrix{Float64})
    cs.fdp .= 0.0
    cs.fx .= fx .+ 0.0
    M = size(fx)[2]
    
    for i in 1:M
        _dodp!(cs.x, cs.fx[:,i], cs.fdp[:,i], cs.npts)
    end

    cs.dirty = false
end



function interp_level(cs::MultiCubicSpline{Float64},
                      x::Float64, i::Int)
    y = zeros(1)
    yp = zeros(1)
    ydp = zeros(1)

    _interp!(x, cs.x, cs.fx[:,i], cs.fdp[:,i], cs.npts,
             1, 0, 0, y, yp, ydp)
    y[1]
end

function(itp::MultiCubicSpline)(x::Real, i::Int, update::Bool=true)
    if update && itp.dirty
        update!(itp, itp.fx)
    end
    interp_level(itp, x, i)
end


