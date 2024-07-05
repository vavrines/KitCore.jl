"""
$(SIGNATURES)

Read text into dictionary

## Arguments
* `filename`: configuration text file
* `allowed`: list of allowed attributes

## Outputs
* dictionary with values of variables
"""
function read_dict(filename::T, allowed) where {T<:AbstractString}
    @info "reading config from $filename"
    println("--------------------------------------------------------------")
    f = open(filename)
    vars = Dict{String,Any}()

    for line in eachline(f)
        # skip comments
        if length(line) == 0 || line[1] == '#'
            continue
        end

        #print("\t");
        #println(line)
        var, val = split(line, "=")
        stripped = strip(var)
        if stripped in allowed
            println(line)

            #vars[stripped] = parse(Float64, val)
            #vars[stripped] = strip(val)
            tmp = tryparse(Float64, val)
            if isa(tmp, Nothing)
                tmp1 = tryparse(Bool, val)
                if isa(tmp1, Nothing)
                    vars[stripped] = strip(val)
                else
                    vars[stripped] = tmp1
                end
            else
                vars[stripped] = isinteger(tmp) ? Int(tmp) : tmp
            end
        end
    end
    println("--------------------------------------------------------------")
    println("")

    return vars
end

"""
$(SIGNATURES)
"""
function read_dict(filename::T) where {T<:AbstractString}
    @info "reading config from $filename"
    println("--------------------------------------------------------------")
    f = open(filename)
    vars = Dict{Symbol,Any}()

    for line in eachline(f)
        if length(line) == 0 || line[1] == '#'
            continue
        end
        println(line)

        var, val = split(line, "=")
        stripped = strip(var)
        stripped = Symbol(stripped)
        val = split(val, "#")[1]

        tmp = tryparse(Float64, val)
        if isa(tmp, Nothing)
            tmp1 = tryparse(Bool, val)
            if isa(tmp1, Nothing)
                vars[stripped] = strip(val)
            else
                vars[stripped] = tmp1
            end
        else
            vars[stripped] = isinteger(tmp) ? Int(tmp) : tmp
        end
    end
    println("--------------------------------------------------------------")
    println("")

    return vars
end
