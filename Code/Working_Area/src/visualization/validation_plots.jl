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





function plot_and_save_validation_full(dict_A, dict_B, separation_dist, N, steps, save_dir)
    println("📊 Generando reporte detallado...")

    # Crear carpetas para organizar
    path_indiv = joinpath(save_dir, "Individual_Z")
    path_pairs = joinpath(save_dir, "Correlations_ZZ")
    mkpath(path_indiv); mkpath(path_pairs)

    all_labels = collect(keys(dict_A))
    z_labels = sort([l for l in all_labels if count(c -> c == 'Z', l) == 1])
    zz_labels = sort([l for l in all_labels if count(c -> c == 'Z', l) == 2])

    # 1. Gráficas individuales ESP para cada qubit
    for lbl in z_labels
        idx = findfirst('Z', lbl)
        p = plot(dict_A[lbl], label="Tray. A", c=:blue, lw=1.5, title="ESP: Qubit $idx")
        plot!(p, dict_B[lbl], label="Tray. B", c=:red, ls=:dash, lw=1.5)
        ylims!(p, (-0.015, 0.015)) 
        savefig(p, joinpath(path_indiv, "ESP_Z$idx.png"))
    end

    # 2. Gráficas para pares ZZ
    for lbl in zz_labels
        indices = [i for (i, c) in enumerate(lbl) if c == 'Z']
        tag = "Z$(indices[1])Z$(indices[2])"
        p = plot(dict_A[lbl], label="Tray. A", c=:blue, lw=1.5, title="ESP Pair: $tag")
        plot!(p, dict_B[lbl], label="Tray. B", c=:red, ls=:dash, lw=1.5)
        #ylims!(p, (-0.015, 0.015)) 
        savefig(p, joinpath(path_pairs, "ESP_$tag.png"))
    end

    # 3. Resumen combinado
    p1 = plot(dict_A[zz_labels[1]], label="Tray. A", c=:blue, title="ESP Sample (ZZ)")
    plot!(p1, dict_B[zz_labels[1]], label="Tray. B", c=:red, ls=:dash)
    
    p2 = bar(separation_dist, label="Separation", c=:green, alpha=0.5)
    
    # Heatmap (Solo sitios individuales para legibilidad)
    h_mat = zeros(N, steps)
    for (i, lbl) in enumerate(z_labels); h_mat[i, :] = dict_A[lbl]; end
    p3 = heatmap(h_mat, c=:viridis, title="Reservoir Map", xlabel="Step", ylabel="Qubit")

    final = plot(p1, p2, p3, layout=@layout([a; b; c]), size=(850, 1000))
    savefig(final, joinpath(save_dir, "Summary_Full_Validation.png"))
end
