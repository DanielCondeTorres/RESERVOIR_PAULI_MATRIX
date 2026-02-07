function get_basis_char(dict_data)
    # Busca qué letra (X, Y, Z) predomina en las llaves
    all_keys = collect(keys(dict_data))
    for char in ['X', 'Y', 'Z']
        if any(k -> contains(k, string(char)), all_keys)
            return char
        end
    end
    return 'Z' # Default
end

"""
1. QUICK VALIDATION PER QUBIT (Auto-Detect)
Genera: Quick_Validation_ESP_Qubit_X1.png, etc.
"""
function plot_quick_validation_per_qubit(separation_dist, dict_A, inputs, N, save_dir)
    basis_char = get_basis_char(dict_A)
    b_str = string(basis_char)
    println("📊 Generando Quick Validation ($b_str) para $N Qubits...")
    
    # Máscaras para scatter
    len = length(inputs); transient = 20
    mask0 = (inputs .== 0) .& (1:len .> transient)
    mask1 = (inputs .== 1) .& (1:len .> transient)
    
    p1 = plot(separation_dist, label="ESP Dist", yscale=:log10, color=:red, 
              title="Global ESP Convergence", xlabel="Step", ylabel="Dist")
    
    for i in 1:N
        # Construye etiqueta dinámica: ej "1X1111"
        lbl_vec = ["1" for _ in 1:N]; lbl_vec[i] = b_str; lbl = join(lbl_vec)
        
        if haskey(dict_A, lbl)
            traj = dict_A[lbl]
            p2 = scatter(findall(mask0), traj[mask0], label="In 0", color=:blue, ms=3, alpha=0.6)
            scatter!(p2, findall(mask1), traj[mask1], label="In 1", color=:orange, ms=3, alpha=0.6,
                     title="Qubit $i (<$lbl>)", ylabel="Exp <$b_str$i>")
            
            p_final = plot(p1, p2, layout=(2,1), size=(800, 800))
            savefig(p_final, joinpath(save_dir, "Quick_Validation_ESP_Qubit_$(b_str)$(i).png"))
        end
    end
end

"""
2. SCATTER MAP ALL QUBITS (Auto-Detect)
Genera: Separability_All_Qubits_Scatter_X.png
"""
function plot_all_qubits_scatter(dict_A, inputs, N, save_dir)    
    basis_char = get_basis_char(dict_A)
    b_str = string(basis_char)
    
    rows = Int(ceil(N / 2))
    p_global = plot(layout=(rows, 2), size=(900, 300*rows), 
                    title="Separability Map - Base $b_str")

    mask0 = (inputs .== 0) .& (1:length(inputs) .> 20)
    mask1 = (inputs .== 1) .& (1:length(inputs) .> 20)

    for i in 1:N
        lbl_vec = ["1" for _ in 1:N]; lbl_vec[i] = b_str; lbl = join(lbl_vec)
        
        if haskey(dict_A, lbl)
            traj = dict_A[lbl]
            scatter!(p_global[i], findall(mask0), traj[mask0], label="In 0", c=:blue, ms=2, alpha=0.6)
            scatter!(p_global[i], findall(mask1), traj[mask1], label="In 1", c=:orange, ms=2, alpha=0.6,
                     title="Q$i ($lbl)", legend=false)
        end
    end
    savefig(p_global, joinpath(save_dir, "Separability_All_Qubits_Scatter_$(b_str).png"))
end

"""
3. BOXPLOTS (Auto-Detect + Sorting Fix)
Genera: Separability_Boxplots_X.png
"""
function plot_separability_boxplots(dict_A, inputs, N, save_dir)
    basis_char = get_basis_char(dict_A)
    b_str = string(basis_char)
    println("📦 Generando Boxplots ($b_str)...")

    # Filtramos solo las etiquetas individuales de la base correcta
    # Y ordenamos por la posición del carácter (para que Q1 salga primero)
    raw_labels = [l for l in keys(dict_A) if count(c -> c == basis_char, l) == 1]
    sorted_labels = sort(raw_labels, by = x -> findfirst(basis_char, x))

    rows = Int(ceil(N/3))
    p = plot(layout=(rows, 3), size=(1000, 350*rows), title="Boxplots Base $b_str")

    for (i, lbl) in enumerate(sorted_labels)
        # Recuperar el índice real del qubit desde la etiqueta
        q_idx = findfirst(basis_char, lbl)
        
        valid_indices = 21:length(inputs) # Skip transient
        data0 = dict_A[lbl][valid_indices][inputs[valid_indices] .== 0]
        data1 = dict_A[lbl][valid_indices][inputs[valid_indices] .== 1]
        
        boxplot!(p[i], [fill("In 0", length(data0)); fill("In 1", length(data1))], 
                 [data0; data1], title="Q$q_idx ($lbl)", legend=false, c=[:blue :red], alpha=0.6)
    end
    savefig(p, joinpath(save_dir, "Separability_Boxplots_$(b_str).png"))
end