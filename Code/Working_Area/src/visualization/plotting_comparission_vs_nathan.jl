using Plots

"""
    plot_validation_comparison(nathan_data, my_data, target_qubits, n_qubits, err_percent)

Genera y guarda una gráfica en la carpeta 'Outputs' del proyecto.
"""
function plot_validation_comparison(nathan_data::Vector{Float64}, 
                                    my_data::Vector{Float64}, 
                                    target_qubits::Vector{Int}, 
                                    n_qubits::Int, 
                                    err_percent::Float64,
                                    obs_char::Char,
                                    limites::Int)
    
    # 1. Construir la etiqueta dinámica (Ej: "IZZII")
    label_chars = fill('I', n_qubits)
    for q_idx in target_qubits
        if 1 <= q_idx <= n_qubits
            label_chars[q_idx] = obs_char
        end
    end
    str_label = String(label_chars)
    
    # 2. Configurar datos
    limites = min(length(nathan_data), length(my_data))
    t = 1:limites
    
    # 3. Graficar
    p = plot(t, nathan_data[1:len], 
             label = "Nathan ($str_label)", 
             color = :black, 
             lw = 2,
             xlabel = "Time Steps",
             ylabel = "Expectation Value",
             title = "Validation: $str_label (Error: $(round(err_percent, digits=2))%)",
             legend = :topright)
             
    plot!(p, t, my_data[1:len], 
          label = "Replica ($str_label)", 
          color = :red, 
          ls = :dash, 
          lw = 2)
    
    # 4. GUARDAR EN LA CARPETA OUTPUTS
    # ---------------------------------------------------------
    # @__DIR__ es: .../Working_Area/src/visualization
    # Subimos 3 niveles: visualization -> src -> Working_Area -> RAÍZ
    root_path = joinpath(@__DIR__, "..", "..", "..")
    output_folder = joinpath(root_path, "Outputs")
    
    # Crear carpeta si no existe (por seguridad)
    if !isdir(output_folder)
        mkpath(output_folder)
    end
    
    filename = "comparison_$(str_label).png"
    full_path = joinpath(output_folder, filename)
    
    savefig(p, full_path)
    println("   💾 Guardado en: $full_path")
    
    return p
end