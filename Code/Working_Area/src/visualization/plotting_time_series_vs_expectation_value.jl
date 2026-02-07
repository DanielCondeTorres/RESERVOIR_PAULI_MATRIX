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




using Plots

"""
    plot_universal(data, limit_x, x_name, y_name, label_user, filename; kwargs...)

Función maestra para graficar series temporales. Permite elegir SCATTER o LINE.

# Argumentos Obligatorios
- `data`: Vector de datos.
- `limit_x`: Límite del eje X.
- `x_name`, `y_name`: Etiquetas de ejes.
- `label_user`: Etiqueta leyenda (si use_qubit_logic=false).
- `filename`: Nombre del archivo de salida.

# Argumentos Opcionales (Keywords)
- `style`: Estilo del gráfico. Opciones: :scatter, :line, :steppost. (Default: :line)
- `use_qubit_logic`: Si true, genera etiqueta tipo "IZXII".
- `target_qubits`, `n_qubits`, `obs_str`: Datos para la etiqueta automática.
"""
function plot_universal(data::Vector{Float64}, 
                        limit_x::Int, 
                        x_name::String, 
                        y_name::String, 
                        label_user::String, 
                        filename::String;
                        # --- Opciones Extra ---
                        style::Symbol = :line,  # <--- NUEVO: :scatter, :line, :steppost
                        use_qubit_logic::Bool = false,
                        target_qubits::Vector{Int} = Int[], 
                        n_qubits::Int = 0, 
                        obs_str::String = "")

    # 1. Lógica de Etiqueta (Flag)
    final_label = label_user
    if use_qubit_logic
        if n_qubits > 0 && obs_str != ""
            char_to_use = obs_str[1] 
            label_chars = fill('I', n_qubits)
            for q_idx in target_qubits
                if 1 <= q_idx <= n_qubits; label_chars[q_idx] = char_to_use; end
            end
            final_label = "Sim ($String(label_chars))"
        end
    end

    # 2. Configuración Estética según el estilo
    # Si es scatter, queremos puntos y no líneas. Si es línea, al revés.
    my_markersize = (style == :scatter) ? 3 : 0
    my_linewidth  = (style == :scatter) ? 0 : 2
    
    # 3. Recorte de datos
    len_real = min(length(data), limit_x)
    t = 1:len_real
    data_cut = data[1:len_real]

    # 4. Graficar
    println("📊 Generando gráfico ($style): $filename ...")
    
    p = plot(t, data_cut, 
             seriestype = style,       # <--- AQUÍ SE APLICA TU ELECCIÓN
             markersize = my_markersize,
             markercolor = :blue,
             markerstrokewidth = 0,
             lw = my_linewidth,        # Grosor de línea
             label = final_label, 
             xlabel = x_name, 
             ylabel = y_name,
             title = "Data: $final_label",
             grid = true,
             legend = :topright,
             # Un poco de margen vertical automático
             ylim = (minimum(data_cut)-0.1, maximum(data_cut)+0.1) 
    )

    # 5. Guardar
    root_path = joinpath(@__DIR__, "..", "..", "..")
    output_folder = joinpath(root_path, "Outputs")
    if !isdir(output_folder); mkpath(output_folder); end
    
    full_path = joinpath(output_folder, filename)
    savefig(p, full_path)
    println("   💾 Guardado en: $full_path")
    
    return p
end