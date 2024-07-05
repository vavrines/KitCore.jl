"""
$(SIGNATURES)

Generate uniform mesh
"""
function uniform_mesh(x0, xnum::Integer, dx)
    points = zeros(xnum)
    for i in eachindex(points)
        points[i] = x0 + (i - 0.5) * dx
    end

    return points
end


"""
$(SIGNATURES)

Equivalent N-dimensional mesh generator as matlab

The difference between `ndgrid` and `meshgrid` is that `meshgrid` produces data that is oriented in Cartesian coordinates (row-major), 
while `ndgrid` produces data related to dimension order (column-major).
"""
ndgrid(v::AV) = copy(v)

"""
$(SIGNATURES)
"""
function ndgrid(v1::AV, v2::AV)
    m, n = length(v1), length(v2)
    v1 = reshape(v1, m, 1)
    v2 = reshape(v2, 1, n)

    return (repeat(v1, 1, n), repeat(v2, m, 1))
end

"""
$(SIGNATURES)
"""
function ndgrid(vs::AV{T}...) where {T}
    ndgrid_fill(a, v, s, snext) = begin
        for j = 1:length(a)
            a[j] = v[div(rem(j - 1, snext), s)+1]
        end
    end

    n = length(vs)
    sz = map(length, vs)
    out = ntuple(i -> Array{T}(undef, sz), n)
    s = 1
    for i = 1:n
        a = out[i]::Array
        v = vs[i]
        snext = s * size(a, i)
        ndgrid_fill(a, v, s, snext)
        s = snext
    end

    return out
end


"""
$(SIGNATURES)

Equivalent structured mesh generator as matlab

The difference between `ndgrid` and `meshgrid` is that `meshgrid` produces data that is oriented in Cartesian coordinates (row-major), 
while `ndgrid` produces data related to dimension order (column-major).
"""
meshgrid(v::AV) = meshgrid(v, v)

"""
$(SIGNATURES)
"""
function meshgrid(x::AV, y::AV)
    X = [i for j in y, i in x]
    Y = [j for j in y, i in x]

    return X, Y
end

"""
$(SIGNATURES)
"""
function meshgrid(x::AV, y::AV, z::AV)
    X = [i for k in z, j in y, i in x]
    Y = [j for k in z, j in y, i in x]
    Z = [k for k in z, j in y, i in x]

    return X, Y, Z
end


"""
$(SIGNATURES)

Find the location index of a point in mesh
"""
function find_idx(x::AV, p; mode = :nonuniform::Symbol)
    if mode == :uniform
        dx = x[2] - x[1]
        return Int(ceil((p - x[1] + 0.5 * dx) / dx)) # point location
    else
        return argmin(abs.(x .- p)) # center location
    end
end


"""
$(SIGNATURES)

Calculate unit normal vector
"""
function unit_normal(p1::T, p2::T) where {T<:AV}
    Δ = p2 .- p1
    l = norm(Δ) + 1e-6

    return [-Δ[2], Δ[1]] ./ l
end

"""
$(SIGNATURES)
"""
function unit_normal(p1::T, p2::T, p3::T) where {T<:AV}
    v1 = p2 .- p1
    v2 = p3 .- p1

    n = cross(v1, v2)
    l = norm(n) + 1e-6

    return n ./ l
end


"""
$(SIGNATURES)

Calculate point-point/line/surface distance
"""
point_distance(p1::T, p2::T) where {T<:AV} = norm(p1 .- p2)

"""
$(SIGNATURES)
"""
function point_distance(p::T, p1::T, p2::T) where {T<:AV}
    x0, y0 = p
    x1, y1 = p1
    x2, y2 = p2

    return abs((x2 - x1) * (y1 - y0) - (x1 - x0) * (y2 - y1)) /
           sqrt((x2 - x1)^2 + (y2 - y1)^2 + 1e-6)
end
