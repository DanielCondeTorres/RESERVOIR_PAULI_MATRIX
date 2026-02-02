using Plots
using Plots
"""
    plot_zz_correlations(steps_vec, expect_dict, folder_path, inputs_vec=nothing)

Busca todas las correlaciones ZZ en el diccionario y genera una gráfica 
en la carpeta del experimento.
"""
function plot_zz_correlations(steps_vec, expect_dict, folder_path, inputs_vec=nothing)
    println("📊 Generando gráfica de correlaciones ZZ...")

    # 1. Filtrar etiquetas con exactamente dos 'Z' (correlaciones de pares)
    zz_labels = sort([k for k in keys(expect_dict) if count(c -> c == 'Z', k) == 2])
    
    if isempty(zz_labels)
        println("⚠️ No se encontraron correlaciones ZZ en el diccionario.")
        return
    end

    # 2. Extraer datos a una matriz
    n_steps = length(steps_vec)
    history_zz = zeros(Float64, n_steps, length(zz_labels))
    for (idx, label) in enumerate(zz_labels)
        history_zz[:, idx] = expect_dict[label]
    end

    # 3. Limpiar etiquetas para la leyenda (ej: "1Z1Z11" -> "Z2Z4")
    # Calculamos la posición real de las Z para que Nathan lo entienda fácil
    clean_labels = String[]
    for l in zz_labels
        indices = [i for (i, char) in enumerate(l) if char == 'Z']
        push!(clean_labels, "Z$(indices[1])Z$(indices[2])")
    end
    clean_labels_reshape = reshape(clean_labels, 1, :)

    # 4. Crear el Plot
    p = plot(steps_vec, history_zz, 
        label = clean_labels_reshape, 
        title = "ZZ Correlations: " * basename(folder_path),
        xlabel = "Step (k)", 
        ylabel = "Expectation <Z_i Z_j>",
        lw = 1.2, 
        palette = :tab20, 
        legend = :outertopright, 
        ylim = (-1.1, 1.1),
        grid = true
    )
    
    # 5. Añadir Input de fondo
    if !isnothing(inputs_vec)
        plot!(p, steps_vec, inputs_vec, label="Input", color=:black, ls=:dash, alpha=0.1, seriestype=:steppre)
    end

    # 6. Guardado
    output_path = joinpath(folder_path, "zz_correlations.png")
    savefig(p, output_path)
    println("✅ Gráfica de correlaciones guardada en: $output_path")
end
function plot_expectation_evolution_easy(steps_vec, results_data, n_sites, folder_path, inputs_vec=nothing)
    println("Generando gráfico de Evolución de Expectación...")

    labels = reshape(["Z$i" for i in 1:n_sites], 1, n_sites)

    p = plot(steps_vec, results_data,
        label = labels,
        xlabel = "Step (k)",
        ylabel = "Expectation <Z_j>",
        title = "Evolution: " * basename(folder_path), # Título dinámico
        lw = 1.5,
        palette = :tab10,
        legend = :outertopright,
        grid = true,
        ylim = (-1.1, 1.1) 
    )
    
    if !isnothing(inputs_vec)
        plot!(p, steps_vec, inputs_vec, label="Input", color=:black, ls=:dash, alpha=0.2, seriestype=:steppre)
    end

    # --- CAMBIO CLAVE AQUÍ ---
    # Usamos folder_path en lugar de SCRIPT_DIR
    output_path = joinpath(folder_path, "expectation_evolution.png")
    savefig(p, output_path)
    println("✅ Gráfico guardado exitosamente en: $output_path")
end

"""
    plot_expectation_evolution(steps_vec, results_data, n_sites)

Genera y guarda el gráfico de la evolución de la expectación <X_j>.
- steps_vec: Eje X (pasos k).
- results_data: Matriz [Pasos x Sitios].
- n_sites: Número de qubits/sitios.
Guarda la imagen en Code/Outputs/expectation_evolution.png
"""
function plot_expectation_evolution(steps_vec::AbstractVector, results_data::Matrix{Float64}, n_sites::Int)
    println("Generando gráfico de Evolución de Expectación...")

    # 1. Etiquetas para la leyenda (Site 1, Site 2...)
    # Usamos reshape para que Plots entienda que son series distintas
    labels = reshape(["Site $(i+1)" for i in 0:n_sites-1], 1, n_sites)

    # 2. Configurar el objeto Plot
    p = plot(steps_vec, results_data,
        label = labels,              # Nombres de las líneas
        xlabel = "Step (k)",
        ylabel = "Expectation <X_j>",
        title = "Evolution of Expectation Value",
        
        # Estilo visual similar al paper
        lw = 2,                      # Grosor de línea
        palette = :tab10,            # Colores categóricos distintos
        legend = :outertopright,     # Leyenda fuera
        grid = true,
        margin = 5Plots.mm,          # Margen extra
        ylim = (-0.1, 1.05)          # Límites verticales para ver bien desde 0 a 1
    )
    
    # Línea de referencia en 0
    hline!(p, [0.0], color=:black, lw=0.8, label="", linestyle=:dot)

    # 3. Guardar el archivo
    # Ruta relativa robusta: sube 3 niveles desde src/visualization/
    output_dir = joinpath(@__DIR__, "..", "..", "..", "Outputs")
    output_dir = abspath(output_dir) # Normaliza la ruta
    
    # Crear carpeta si no existe
    if !isdir(output_dir)
        mkpath(output_dir)
    end
    
    file_path = joinpath(output_dir, "expectation_evolution.png")
    savefig(p, file_path)
    
    println("✅ Gráfico guardado exitosamente en: $file_path")
    # display(p) # Descomenta si quieres verlo en pantalla al ejecutar
end




