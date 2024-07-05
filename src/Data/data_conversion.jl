"""
$(SIGNATURES)

Extract subarray except the last column
"""
function extract_last(a::AA{T}, idx::Integer; mode = :view::Symbol) where {T}
    if mode == :copy

        if ndims(a) == 2
            sw = a[:, idx]
        elseif ndims(a) == 3
            sw = a[:, :, idx]
        elseif ndims(a) == 4
            sw = a[:, :, :, idx]
        elseif ndims(a) == 5
            sw = a[:, :, :, :, idx]
        end

    elseif mode == :view

        if ndims(a) == 2
            sw = @view a[:, idx]
        elseif ndims(a) == 3
            sw = @view a[:, :, idx]
        elseif ndims(a) == 4
            sw = @view a[:, :, :, idx]
        elseif ndims(a) == 5
            sw = @view a[:, :, :, :, idx]
        end

    else

        throw("Error in extraction mode setup")

    end

    return sw
end
