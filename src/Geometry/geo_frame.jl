"""
$(SIGNATURES)

Transform local flow variables to global frame
"""
function global_frame(w::AV, cosa, sina)
    if eltype(w) <: Int
        G = similar(w, Float64)
    else
        G = similar(w)
    end

    if length(w) == 2
        G[1] = w[1] * cosa - w[2] * sina
        G[2] = w[1] * sina + w[2] * cosa
    elseif length(w) == 4
        G[1] = w[1]
        G[2] = w[2] * cosa - w[3] * sina
        G[3] = w[2] * sina + w[3] * cosa
        G[4] = w[4]
    else
        throw("local -> global: dimension dismatch")
    end

    return G
end

global_frame(w::AV, n::AV) = global_frame(w, n[1], n[2])

"""
$(SIGNATURES)
"""
function global_frame(w::AV, dirccos::AM)
    if eltype(w) <: Int
        G = similar(w, Float64)
    else
        G = similar(w)
    end

    if length(w) == 3
        G[1] = w[1] * dirccos[1, 1] + w[2] * dirccos[2, 1] + w[3] * dirccos[3, 1]
        G[2] = w[1] * dirccos[1, 2] + w[2] * dirccos[2, 2] + w[3] * dirccos[3, 2]
        G[3] = w[1] * dirccos[1, 3] + w[2] * dirccos[2, 3] + w[3] * dirccos[3, 3]
    elseif length(w) == 5
        G[1] = w[1]
        G[2] = w[2] * dirccos[1, 1] + w[3] * dirccos[2, 1] + w[4] * dirccos[3, 1]
        G[3] = w[2] * dirccos[1, 2] + w[3] * dirccos[2, 2] + w[4] * dirccos[3, 2]
        G[4] = w[2] * dirccos[1, 3] + w[3] * dirccos[2, 3] + w[4] * dirccos[3, 3]
        G[5] = w[5]
    else
        throw("local -> global: dimension dismatch")
    end

    return G
end


"""
$(SIGNATURES)

Transform global flow variables to local frame
"""
function local_frame(w::AV, cosa, sina)
    if eltype(w) <: Int
        L = similar(w, Float64)
    else
        L = similar(w)
    end

    if length(w) == 2
        L[1] = w[1] * cosa + w[2] * sina
        L[2] = w[2] * cosa - w[1] * sina
    elseif length(w) == 4
        L[1] = w[1]
        L[2] = w[2] * cosa + w[3] * sina
        L[3] = w[3] * cosa - w[2] * sina
        L[4] = w[4]
    else
        throw("global -> local: dimension dismatch")
    end

    return L
end

"""
$(SIGNATURES)
"""
local_frame(w::AV, n) = local_frame(w::AV, n[1], n[2])

"""
$(SIGNATURES)
"""
function local_frame(w::AV, dirccos::AM)
    if eltype(w) <: Int
        L = similar(w, Float64)
    else
        L = similar(w)
    end

    if length(w) == 3
        L[1] = w[1] * dirccos[1, 1] + w[2] * dirccos[1, 2] + w[3] * dirccos[1, 3]
        L[2] = w[1] * dirccos[2, 1] + w[2] * dirccos[2, 2] + w[3] * dirccos[2, 3]
        L[3] = w[1] * dirccos[3, 1] + w[2] * dirccos[3, 2] + w[3] * dirccos[3, 3]
    elseif length(w) == 5
        L[1] = w[1]
        L[2] = w[2] * dirccos[1, 1] + w[3] * dirccos[1, 2] + w[4] * dirccos[1, 3]
        L[3] = w[2] * dirccos[2, 1] + w[3] * dirccos[2, 2] + w[4] * dirccos[2, 3]
        L[4] = w[2] * dirccos[3, 1] + w[3] * dirccos[3, 2] + w[4] * dirccos[3, 3]
        L[5] = w[5]
    else
        throw("global -> local: dimension dismatch")
    end

    return L
end


"""
$(SIGNATURES)

Transform velocity to local frame
"""
function local_velocity(u, v, cosa, sina)
    vn = @. u * cosa + v * sina
    vt = @. v * cosa - u * sina

    return vn, vt
end

local_velocity(u, v, n) = local_velocity(u, v, n[1], n[2])


"""
$(SIGNATURES)

Transform velocity to global frame
"""
function global_velocity(u, v, cosa, sina)
    vn = @. u * cosa - v * sina
    vt = @. u * sina + v * cosa

    return vn, vt
end

global_velocity(u, v, n) = global_velocity(u, v, n[1], n[2])
