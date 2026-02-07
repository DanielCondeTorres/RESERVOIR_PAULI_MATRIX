"""
    plot_quick_validation_per_qubit(separation_dist, dict_A, inputs, N, save_dir)

Genera una imagen por cada Qubit (del 1 al N) mostrando:
1. (Arriba) La convergencia global ESP (Distancia entre trayectorias).
2. (Abajo) La separabilidad específica de ESE qubit (Scatter plot).
"""
function plot_quick_validation_per_qubit(separation_dist, dict_A, inputs, N, save_dir)
    println("📊 Generando Quick Validation detallada para cada Qubit (1 a $N)...")
    
    # 1. Pre-calcular máscaras (Ignoramos los primeros 20 pasos de transitorio)
    # Esto se hace una sola vez para ganar velocidad
    len = length(inputs)
    transient = 20
    mask0 = (inputs .== 0) .& (1:len .> transient)
    mask1 = (inputs .== 1) .& (1:len .> transient)
    
    # 2. Plot 1: ESP Distance (Es idéntico para todos los qubits)
    # Lo definimos fuera del bucle para reutilizarlo
    p1 = plot(separation_dist, label="ESP Distance", yscale=:log10, 
              color=:red, lw=2, title="ESP Convergence (Global)",
              xlabel="Steps", ylabel="||Rho_A - Rho_B||", 
              legend=:topright, grid=true, gridalpha=0.3)
    
    # 3. Bucle para generar Scatter de cada Qubit individual
    for i in 1:N
        # Construimos la etiqueta dinámica: ej. i=2, N=4 -> "1Z11"
        lbl_vec = ["1" for _ in 1:N]
        lbl_vec[i] = "Z"
        lbl = join(lbl_vec)
        
        if haskey(dict_A, lbl)
            traj = dict_A[lbl]
            
            # Plot 2: Scatter específico del Qubit i
            p2 = scatter(findall(mask0), traj[mask0], label="Input 0 (|0>)", 
                         color=:blue, ms=3, alpha=0.7, markerstrokewidth=0)
            scatter!(p2, findall(mask1), traj[mask1], label="Input 1 (|+>)", 
                     color=:orange, ms=3, alpha=0.7, markerstrokewidth=0,
                     title="Separabilidad Qubit $i (<$lbl>)", 
                     xlabel="Step", ylabel="Expectation <Z$i>")
            
            # Combinar ESP + Scatter
            p_final = plot(p1, p2, layout=(2,1), size=(800, 800))
            
            # Guardar archivo único para este qubit
            filename = "Quick_Validation_ESP_Qubit_$(i).png"
            output_path = joinpath(save_dir, filename)
            savefig(p_final, output_path)
            
            # Solo hacemos display del primero para no saturar la pantalla
            if i == 1; display(p_final); end
        else
            println("⚠️ Aviso: No se encontraron datos para el qubit $i (Etiqueta: $lbl)")
        end
    end
    
    println("✅ ¡Imágenes generadas! Revisa la carpeta: $save_dir")
end


function plot_all_qubits_scatter(dict_A, inputs, N, save_dir)    
    # Calculamos layout (ej: 2 columnas)
    cols = 2
    rows = Int(ceil(N / cols))
    p_global = plot(layout=(rows, cols), size=(900, 300*rows), title="Separability Map (All Qubits)")

    for i in 1:N
        # Construir etiqueta: "111Z11" para el qubit i
        lbl_vec = ["1" for _ in 1:N]
        lbl_vec[i] = "Z"
        lbl = join(lbl_vec)
        
        if haskey(dict_A, lbl)
            traj = dict_A[lbl]
            # Filtramos el transitorio inicial (> 20) para ver la separabilidad real
            mask0 = (inputs .== 0) .& (1:length(inputs) .> 20)
            mask1 = (inputs .== 1) .& (1:length(inputs) .> 20)
            
            scatter!(p_global[i], findall(mask0), traj[mask0], label="In 0", color=:blue, ms=2, alpha=0.6)
            scatter!(p_global[i], findall(mask1), traj[mask1], label="In 1", color=:orange, ms=2, alpha=0.6,
                     title="Qubit $i", ylabel="<$lbl>", legend=(i==1))
        end
    end
    
    savefig(p_global, joinpath(save_dir, "Separability_All_Qubits_Scatter.png"))
    display(p_global)
end