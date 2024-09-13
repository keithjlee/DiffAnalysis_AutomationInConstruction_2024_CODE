function find_bands(bool_vector::Vector{Bool})
    true_band_starts = []
    true_band_ends = []
    false_band_starts = []
    false_band_ends = []
    
    current_band_start = 1
    current_band_value = bool_vector[1]
    
    for i in 2:length(bool_vector)
        if bool_vector[i] != current_band_value
            if current_band_value
                push!(true_band_starts, current_band_start)
                push!(true_band_ends, i - 1)
            else
                push!(false_band_starts, current_band_start)
                push!(false_band_ends, i - 1)
            end
            current_band_start = i
            current_band_value = bool_vector[i]
        end
    end
    
    # Handle the last band
    if current_band_value
        push!(true_band_starts, current_band_start)
        push!(true_band_ends, length(bool_vector))
    else
        push!(false_band_starts, current_band_start)
        push!(false_band_ends, length(bool_vector))
    end
    
    return (true_band_starts, true_band_ends, false_band_starts, false_band_ends)
end

function check_cstr(cstr, baseline, tol = 1e-6)
    all(cstr ./ baseline .< tol)
end