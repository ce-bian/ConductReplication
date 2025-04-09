function texout_new(data::Any, filepath::String; rounding::Union{Nothing, Vector{Int}}=nothing, col_select::Union{Nothing, Vector{Symbol}, Vector{Int}} = nothing)
    if isa(data, Matrix)
        data = DataFrame(data, :auto)
    end
    if !isnothing(col_select)
        # Select the desired columns if specified
        data = select(data, col_select)
    end
    if isnothing(rounding)
        # Default to 3 digits if not specified
        rounding = fill(3, ncol(data))
    end
    # Ensure the rounding vector matches the number of columns
    if length(rounding) != ncol(data)
        error("The length of the rounding vector must match the number of columns in the data.")
    end
    open(filepath, "w") do io
        for row in eachrow(data)
            # Apply specific rounding for each column based on the rounding vector
            formatted_row = [
                if isa(row[i], AbstractFloat)
                    rounding[i] == 0 ? round(Int, row[i]) : round(row[i], digits=rounding[i])
                else
                    row[i]
                end
                for i in 1:ncol(data)
            ]
            println(io, join(formatted_row, " & "), " \\\\")
        end
    end
end