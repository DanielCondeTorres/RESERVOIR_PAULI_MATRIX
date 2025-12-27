using Plots

function plot_validation_comparison(nathan_data::Vector{Float64}, 
    my_data::Vector{Float64}, 
    target_qubits::Vector{Int}, 
    n_qubits::Int, 
    err_percent::Float64,
    obs_str::String,      # <--- AHORA ACEPTA STRING "Z"
    limites::Int)

    # 1. Construir la etiqueta dinámica
    # Extraemos el primer carácter del string (ej: de "Z" sacamos 'Z')
    char_to_use = obs_str[1] 
    
    label_chars = fill('I', n_qubits)
    for q_idx in target_qubits
        if 1 <= q_idx <= n_qubits
            label_chars[q_idx] = char_to_use
        end
    end
    str_label = String(label_chars)
    
    # 2. Configurar datos
    len_real = min(length(nathan_data), length(my_data), limites)
    t = 1:len_real
    
    # 3. Graficar
    p = plot(t, nathan_data[1:len_real], 
             label = "Nathan ($str_label)", 
             color = :black, 
             lw = 2,
             xlabel = "Time Steps",
             ylabel = "Expectation Value",
             title = "Validation: $str_label (Error: $(round(err_percent, digits=2))%)",
             legend = :topright)
             
    plot!(p, t, my_data[1:len_real], 
          label = "Replica ($str_label)", 
          color = :red, 
          ls = :dash, 
          lw = 2)
    
    # 4. Guardar
    root_path = joinpath(@__DIR__, "..", "..", "..")
    output_folder = joinpath(root_path, "Outputs")
    
    if !isdir(output_folder); mkpath(output_folder); end
    
    filename = "comparison_$(str_label).png"
    full_path = joinpath(output_folder, filename)
    
    savefig(p, full_path)
    println("   💾 Guardado en: $full_path")
    
    return p
end