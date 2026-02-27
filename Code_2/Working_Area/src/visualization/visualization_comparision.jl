function plot_and_save_validation_full(dict_A, dict_B,  N, steps, save_dir; label_A="Tray. A", label_B="Tray. B")
    println("📊 Generando reporte detallado ($label_A vs $label_B)...")
    all_labels = collect(keys(dict_A)) 
    # 1. Detectar qué base estamos usando (X, Y o Z)
    basis_char = 'Z' # Default
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
    path_triads = joinpath(save_dir, "Correlations_$(b_str)$(b_str)$(b_str)") # <--- NUEVA CARPETA
    mkpath(path_indiv); mkpath(path_pairs); mkpath(path_triads)
    # 3. Filtrar etiquetas
    # Individuales (1 letra), Pares (2 letras), Tríadas (3 letras)
    z_labels = sort([l for l in all_labels if count(c -> c == basis_char, l) == 1])
    zz_labels = sort([l for l in all_labels if count(c -> c == basis_char, l) == 2])
    zzz_labels = sort([l for l in all_labels if count(c -> c == basis_char, l) == 3]) 
    # 4. Gráficas individuales
    for lbl in z_labels
        idx = findfirst(isequal(basis_char), lbl) 
        p = plot(dict_A[lbl], label=label_A, c=:blue, lw=1.5, title="ESP: Qubit $idx ($b_str)")
        plot!(p, dict_B[lbl], label=label_B, c=:red, ls=:dash, lw=1.5)
        ylims!(p, -1.1, 1.1)
        savefig(p, joinpath(path_indiv, "ESP_$(b_str)$idx.png"))
    end

    # 5. Gráficas para pares
    for lbl in zz_labels
        indices = [i for (i, c) in enumerate(lbl) if c == basis_char]
        if length(indices) >= 2
            tag = "$(b_str)$(indices[1])$(b_str)$(indices[2])"
            p = plot(dict_A[lbl], label=label_A, c=:blue, lw=1.5, title="ESP Pair: $tag")
            plot!(p, dict_B[lbl], label=label_B, c=:red, ls=:dash, lw=1.5)
            ylims!(p, -1.1, 1.1)
            savefig(p, joinpath(path_pairs, "ESP_$tag.png"))
        end
    end

    # 6. Gráficas para Tríadas (Xi Xj Xk)
    for lbl in zzz_labels
        indices = [i for (i, c) in enumerate(lbl) if c == basis_char]
        if length(indices) >= 3
            tag = "$(b_str)$(indices[1])$(b_str)$(indices[2])$(b_str)$(indices[3])"
            p = plot(dict_A[lbl], label=label_A, c=:blue, lw=1.0, title="ESP Triad: $tag")
            plot!(p, dict_B[lbl], label=label_B, c=:red, ls=:dash, lw=1.0)
            ylims!(p, -1.1, 1.1)
            savefig(p, joinpath(path_triads, "ESP_$tag.png"))
        end
    end

    # 7. Resumen combinado
    p1 = plot(title="Sample Comparison")
    if !isempty(z_labels)
        sample_lbl = z_labels[1]
        p1 = plot(dict_A[sample_lbl], label=label_A, c=:blue, title="Sample ($sample_lbl)")
        plot!(p1, dict_B[sample_lbl], label=label_B, c=:red, ls=:dash)
        ylims!(p1, -1.1, 1.1)
    end
end