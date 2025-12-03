using Plots

"""
    plot_stm_capacity(capacities, max_delay, n_qubits)

Grafica la Capacidad de Memoria a Corto Plazo (STM) vs el Retraso (tau).
Guarda la imagen en: Code/Outputs/stm_capacity_final.png
"""
function plot_stm_capacity(capacities::Vector{Float64}, max_delay::Int, n_qubits::Int)
    println("Generando gráfico de Capacidad STM...")

    # Configuración del gráfico
    p = plot(1:max_delay, capacities, 
         marker = :o, 
         label = "STM Capacity", 
         xlabel = "Delay (tau)", 
         ylabel = "Capacity C",
         title = "Short Term Memory Capacity (N=$n_qubits)",
         ylim = (0, 1.1),
         lw = 2,
         color = :blue,
         grid = true,
         legend = :topright,
         margin = 5Plots.mm
    )

    # Gestión de rutas para guardar
    # Sube niveles desde src/visualization/ hasta llegar a Code/Outputs
    output_dir = joinpath(@__DIR__, "..", "..", "..", "Outputs")
    output_dir = abspath(output_dir)
    
    if !isdir(output_dir)
        mkpath(output_dir)
    end
    
    file_path = joinpath(output_dir, "stm_capacity_final.png")
    savefig(p, file_path)
    
    println("✅ Gráfico guardado exitosamente en: $file_path")
end