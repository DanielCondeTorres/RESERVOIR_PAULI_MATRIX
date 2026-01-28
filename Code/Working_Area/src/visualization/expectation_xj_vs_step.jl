using Plots
function plot_expectation_evolution_easy(steps_vec, results_data, n_sites, inputs_vec=nothing)
    println("Generando gráfico de Evolución de Expectación...")

    # Etiquetas correctas (Vector fila 1xN)
    labels = reshape(["Z$i" for i in 1:n_sites], 1, n_sites)

    p = plot(steps_vec, results_data,
        label = labels,
        xlabel = "Step (k)",
        ylabel = "Expectation <Z_j>",
        title = "Evolution of Expectation Value",
        lw = 1.5,
        palette = :tab10,
        legend = :outertopright,
        grid = true,
        # Ajustado a (-1.1, 1.1) porque <Z> puede ser negativo
        ylim = (-1.1, 1.1) 
    )
    
    # Si pasamos los inputs, los graficamos en el fondo
    if !isnothing(inputs_vec)
        plot!(p, steps_vec, inputs_vec, label="Input", color=:black, ls=:dash, alpha=0.2, seriestype=:steppre)
    end

    # Guardado robusto
    output_path = joinpath(SCRIPT_DIR, "expectation_evolution.png")
    savefig(p, output_path)
    println("✅ Gráfico guardado en: $output_path")
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




