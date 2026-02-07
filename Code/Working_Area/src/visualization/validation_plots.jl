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
    println("📊 Generando reporte detallado (Corrección de nombres X/Y)...")

    all_labels = collect(keys(dict_A))
    
    # 1. Detectar qué base estamos usando (X, Y o Z)
    # Buscamos en las etiquetas cuál de estas letras aparece
    basis_char = 'Z' 
    for char in ['X', 'Y', 'Z']
        if any(k -> contains(k, string(char)), all_labels)
            basis_char = char
            break
        end
    end
    b_str = string(basis_char)

    # 2. Crear carpetas con nombres dinámicos
    path_indiv = joinpath(save_dir, "Individual_$(b_str)")
    path_pairs = joinpath(save_dir, "Correlations_$(b_str)$(b_str)")
    mkpath(path_indiv); mkpath(path_pairs)

    # 3. Filtrar etiquetas usando la base detectada
    z_labels = sort([l for l in all_labels if count(c -> c == basis_char, l) == 1])
    zz_labels = sort([l for l in all_labels if count(c -> c == basis_char, l) == 2])

    # 4. Gráficas individuales (Ej: X1, X2...)
    for lbl in z_labels
        idx = findfirst(basis_char, lbl)
        # Título y nombre de archivo dinámico
        p = plot(dict_A[lbl], label="Tray. A", c=:blue, lw=1.5, title="ESP: Qubit $idx ($b_str)")
        plot!(p, dict_B[lbl], label="Tray. B", c=:red, ls=:dash, lw=1.5)
        savefig(p, joinpath(path_indiv, "ESP_$(b_str)$idx.png"))
    end

    # 5. Gráficas para pares (Ej: X4X5, X4X6...)
    for lbl in zz_labels
        indices = [i for (i, c) in enumerate(lbl) if c == basis_char]
        # AQUÍ ESTABA EL ERROR: Ahora usa b_str en lugar de "Z"
        tag = "$(b_str)$(indices[1])$(b_str)$(indices[2])"
        
        p = plot(dict_A[lbl], label="Tray. A", c=:blue, lw=1.5, title="ESP Pair: $tag")
        plot!(p, dict_B[lbl], label="Tray. B", c=:red, ls=:dash, lw=1.5)
        savefig(p, joinpath(path_pairs, "ESP_$tag.png"))
    end

    # 6. Resumen combinado
    # Tomamos el primer par de la lista para el ejemplo
    if !isempty(zz_labels)
        p1 = plot(dict_A[zz_labels[1]], label="Tray. A", c=:blue, title="ESP Sample ($(zz_labels[1]))")
        plot!(p1, dict_B[zz_labels[1]], label="Tray. B", c=:red, ls=:dash)
    else
        p1 = plot(title="No pairs found")
    end
    
    p2 = plot(separation_dist, label="Separation", c=:green, lw=2, title="Global ESP Distance")
    
    # Heatmap (Ordenado por qubit)
    h_mat = zeros(N, steps)
    for (i, lbl) in enumerate(z_labels)
        q_idx = findfirst(basis_char, lbl)
        if q_idx <= N; h_mat[q_idx, :] = dict_A[lbl]; end
    end
    p3 = heatmap(h_mat, c=:viridis, title="Reservoir Map ($b_str)", xlabel="Step", ylabel="Qubit")

    final = plot(p1, p2, p3, layout=@layout([a; b; c]), size=(850, 1000))
    savefig(final, joinpath(save_dir, "Summary_Full_Validation_$(b_str).png"))
end