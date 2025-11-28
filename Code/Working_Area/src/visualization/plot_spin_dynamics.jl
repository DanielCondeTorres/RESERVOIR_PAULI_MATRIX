using Plots

"""
    plot_mixed_dynamics(steps, results, labels)
    Grafica la dinámica mixta (Sitios y Enlaces) con etiquetas personalizadas.
    Guarda en: Code/Outputs/spin_dynamics_mixed.png
"""
function plot_mixed_dynamics(steps::Vector{Float64}, results::Matrix{Float64}, labels::Matrix{String})
    println("Generando gráfico de Dinámica Mixta (Sitios + Enlaces)...")

    p = plot(steps, results,
        label = labels,
        xlabel = "Input Step (k)",
        ylabel = "Amplitude",
        title = "Propagación Detallada (Z_i y Enlaces)",
        lw = 2,
        alpha = 0.8,
        palette = :tab10,
        legend = :outertopright,
        grid = true,
        ylim = (-1.1, 1.1), # Rango completo para ver ondas negativas
        size = (900, 500),
        margin = 5Plots.mm
    )
    hline!(p, [0.0], color=:black, lw=1, label="", linestyle=:dot)

    # Guardar
    output_dir = joinpath(@__DIR__, "..", "..", "..", "Outputs")
    output_dir = abspath(output_dir)
    if !isdir(output_dir); mkpath(output_dir); end
    
    file_path = joinpath(output_dir, "spin_dynamics_mixed.png")
    savefig(p, file_path)
    println("✅ Gráfico guardado exitosamente en: $file_path")
end