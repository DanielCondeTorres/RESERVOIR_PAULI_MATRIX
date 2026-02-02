using Plots

"""
    plot_and_save_validation(hist_A, hist_B, separation_dist, N, steps, save_dir)

Genera y guarda:
1. Un resumen gráfico combinado (Heatmap + Separación + ESP Qubit 2).
2. Gráficas individuales de convergencia ESP para cada qubit (Z1, Z2... ZN).
"""
function plot_and_save_validation(hist_A, hist_B, separation_dist, N, steps, save_dir)
    println("📊 Generando gráficas de validación...")

    # ==========================================================================
    # 1. GENERAR GRÁFICAS INDIVIDUALES POR SITE (QUBIT)
    # ==========================================================================
    # Esto te permite ver cómo la propiedad de eco (ESP) varía según la distancia al input
    plots_dir = joinpath(save_dir, "Individual_Qubits")
    mkpath(plots_dir)

    for i in 1:N
        p = plot(hist_A[:, i], label="Tray. A", c=:blue, lw=1.5, title="ESP: Qubit $i")
        plot!(p, hist_B[:, i], label="Tray. B (Rotada)", c=:red, ls=:dash, lw=1.5)
        ylabel!(p, "<Z$i>")
        xlabel!(p, "Step (k)")
        ylims!(p, (-1.1, 1.1)) # Fijar límites para ver bien la escala
        
        savefig(p, joinpath(plots_dir, "ESP_Convergence_Qubit_$i.png"))
    end
    println("   ➜ Gráficas individuales guardadas en: $plots_dir/")

    # ==========================================================================
    # 2. GENERAR RESUMEN COMBINADO (SUMMARY)
    # ==========================================================================
    # Usamos el Qubit 2 para la gráfica de ESP del resumen, ya que suele mostrar
    # mejor la convergencia que el Qubit 1 (que se borra en cada paso).
    
    # A. ESP (Qubit 2)
    p1 = plot(hist_A[:, 2], label="Tray. A (Qubit 2)", c=:blue, lw=1.2)
    plot!(p1, hist_B[:, 2], label="Tray. B (Qubit 2)", c=:red, ls=:dash, lw=1.2)
    title!(p1, "ESP Validation (Qubit 2 - Indirect Memory)")
    ylabel!(p1, "<Z2>")

    # B. Separación
    p2 = bar(separation_dist, label="Distancia (Euclídea)", c=:green, alpha=0.6, lw=0)
    title!(p2, "Separation: Sensitivity to Input Change")
    ylabel!(p2, "Distance")

    # C. Heatmap
    p3 = heatmap(1:steps, 1:N, hist_A', c=:viridis, title="Reservoir Heatmap")
    ylabel!(p3, "Qubit Index")
    xlabel!(p3, "Step (k)")
    
    # Layout Final
    l = @layout [a; b; c]
    final_plot = plot(p1, p2, p3, layout=l, size=(800, 1000), margin=5Plots.mm)
    
    # Guardar Resumen
    summary_path = joinpath(save_dir, "Summary_Validation.png")
    savefig(final_plot, summary_path)
    
    println("   ➜ Resumen guardado en: $summary_path")
end