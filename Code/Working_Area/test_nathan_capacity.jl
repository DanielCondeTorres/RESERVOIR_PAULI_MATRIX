using LinearAlgebra
using Plots
using Random
using JLD2
using FileIO
using Statistics
function build_operator_from_string(op_str::String)
    # Definimos las matrices básicas
    sx = [0. 1.; 1. 0.]; sy = [0. -im; im 0.]; sz = [1. 0.; 0. -1.]; id = [1. 0.; 0. 1.]
    mapping = Dict('X'=>sx, 'Y'=>sy, 'Z'=>sz, 'I'=>id, '1'=>id) # '1' suele ser identidad en Nathan
    
    # Construimos el producto de Kronecker
    char_list = collect(op_str)
    if length(char_list) == 0; return nothing; end
    
    # Primer qubit
    c1 = uppercase(char_list[1])
    op = haskey(mapping, c1) ? mapping[c1] : id
    
    # Resto de qubits
    for c in char_list[2:end]
        c_upper = uppercase(c)
        mat = haskey(mapping, c_upper) ? mapping[c_upper] : id
        op = kron(op, mat)
    end
    return op
end
# ==============================================================================
# 1. CONFIGURACIÓN DE RUTAS Y CARGA DE MÓDULOS
# ==============================================================================
const SCRIPT_DIR = @__DIR__
function include_rel(path...)
    include(joinpath(SCRIPT_DIR, path...))
end

try
    include_rel("src/operator_terms/pauli_algebra.jl")
    include_rel("src/operator_terms/hamiltonian.jl")
    include_rel("src/utils/initial_state.jl")
    include_rel("src/utils/measurements.jl")
    include_rel("src/utils/dynamics.jl")
    include_rel("src/auxiliary_scripts/save_files_like_nathan.jl")
    include_rel("src/utils/injection_EraseWrite.jl")
    include_rel("src/visualization/validation_plots.jl")
    include_rel("src/utils/quantum_channels.jl")
    include_rel("src/metrics/separability_metrics.jl")
    include_rel("src/visualization/separability_plots.jl")
    include_rel("src/visualization/separability_and_Echo_property.jl")
    include_rel("src/capacity_training/qrc_training.jl")
    include_rel("src/capacity_training/calculate_stm_capacity.jl")
catch e
    println("⚠️  Nota: Algunos módulos no cargaron. Asegúrate de que las rutas 'src/...' sean correctas.")
end

# ==============================================================================
# 2. FUNCIONES DE DEPHASING (Versión Corregida)
# ==============================================================================
# ==============================================================================
# 3. PARÁMETROS DEL EXPERIMENTO
# ==============================================================================
N = 6
T_evol = 1.0        
h_val = 1.0
gamma = 0.04      
n_substeps = 100
dt = T_evol / n_substeps
projective_mode = false
gamma_value = 0.8

# Configuración de Archivos
INPUT_FILE = "6_3_2_all_zeros_12345.jld2"
Experiment_name = "Schrodinger_Replay_s_vec"
Experiment_path = joinpath(SCRIPT_DIR, "../$Experiment_name")
if !isdir(Experiment_path); mkpath(Experiment_path); end

# Parámetros STM
WASHOUT = 50      
MAX_DELAY = 15     

# ==============================================================================
# 4. SIMULACIÓN PRINCIPAL
# ==============================================================================
function run_schrodinger_esp_task()
    println("🚀 Iniciando experimento: $Experiment_name")
    
    # --- CARGA DE INPUTS (Buscando s_vec) ---
    input_file_path = isfile(INPUT_FILE) ? INPUT_FILE : joinpath(SCRIPT_DIR, "../Input_Data", INPUT_FILE)
    if !isfile(input_file_path); error("❌ No se encuentra el archivo de inputs: $input_file_path"); end
    
    println("📂 Cargando inputs desde: $input_file_path")
    
    data_loaded = load(input_file_path)
    meta = data_loaded["meta_data_dict"]
    nathan_results = data_loaded["expect_dict"]
    
    # --- SETUP SISTEMA (Usando parámetros del archivo) ---
    N = Int(meta["N"])
    h_file = Float64(meta["h"])
    g_file = Float64(meta["g"]) 
    J_vec_file = Vector{Float64}(meta["Jvec"])
    inputs = Vector{Float64}(meta["s_vec"])
    steps = Int(meta["num_steps"])
    dt = meta["dt"]

    println("   ✅ Parámetros cargados: N=$N, h=$h_file, g=$g_file, Steps=$steps")

    # --- C. CONSTRUCCIÓN ---
    H_op = hamiltonian_nathan_XX(N, J_vec_file, h_file)
    #H_op = build_nathan_all_to_all_XX(N, h_val)
    H_dense = operator_to_dense_matrix(H_op, N)
    rho_A = operator_to_dense_matrix(initial_state_all_zeros(N), N)
    rho_B = operator_to_dense_matrix(initial_state_all_ones(N), N)

    # --- DEFINICIÓN DE OBSERVABLES (HÍBRIDA) ---
    obs_matrices = Dict{String, Matrix{ComplexF64}}()
    
    # 1. Añadimos tus observables estándar (1-body, 2-body)
    all_labels = String[]
    sx=[0. 1.; 1. 0.]; sy=[0. -im; im 0.]; sz=[1. 0.; 0. -1.]; id=[1. 0.; 0. 1.]
    for (b_char, b_op) in zip(['X','Y','Z'], [sx, sy, sz])
        b_str = string(b_char)
        for i in 1:N
            lbl = ["1" for _ in 1:N]; lbl[i] = b_str; s_lbl = join(lbl)
            obs_matrices[s_lbl] = build_operator_from_string(s_lbl)
        end
        for i in 1:N, j in (i+1):N
            lbl = ["1" for _ in 1:N]; lbl[i]=b_str; lbl[j]=b_str; s_lbl = join(lbl)
            obs_matrices[s_lbl] = build_operator_from_string(s_lbl)
        end
    end

    # 2. IMPORTANTE: Añadimos TODO lo que haya en el archivo de Nathan
    #    Esto asegura que si él tiene "Z1Z111", nosotros también lo medimos.
    println("   🔍 Escaneando observables de Nathan...")
    count_new = 0
    for key in keys(nathan_results)
        # Solo procesamos si la clave parece un operador (longitud N)
        if length(key) == N && !haskey(obs_matrices, key)
            try
                op_mat = build_operator_from_string(key)
                obs_matrices[key] = op_mat
                count_new += 1
            catch
                # Si falla (no es un string de Pauli válido), lo ignoramos
            end
        end
    end
    println("   ✅ Se han añadido $count_new observables extra del archivo original.")

    # Inicializamos diccionarios
    # OJO: Usamos keys(obs_matrices) para cubrir tanto los tuyos como los de Nathan
    dict_A = Dict(lbl => zeros(Float64, steps) for lbl in keys(obs_matrices))
    dict_B = Dict(lbl => zeros(Float64, steps) for lbl in keys(obs_matrices))
    separation_dist = zeros(Float64, steps)

    # --- BUCLE DE DINÁMICA ---
    println("🔄 Ejecutando dinámica...")
    for k in 1:steps
        s_k = Float64(inputs[k])
        rz = 1.0 - 2.0 * s_k; rx = 0.0; ry = 0.0 
        
        rho_A = inject_state_EraseWrite_matrix(rho_A, 0, rz, rx, ry)
        rho_B = inject_state_EraseWrite_matrix(rho_B, 0, rz, rx, ry)

        for _ in 1:n_substeps
            rho_A = step_rk4_matrix(rho_A, H_dense, dt)
            rho_B = step_rk4_matrix(rho_B, H_dense, dt)
        end
        
        rho_A = apply_global_dephasing_matrix(rho_A, g_file, "Z")
        rho_B = apply_global_dephasing_matrix(rho_B, g_file, "Z")

        d_sq = 0.0
        for (lbl, op) in obs_matrices
            val_A, rho_A_new = measure_observable_matrix(rho_A, op, projective=projective_mode, gamma_value)
            val_B, rho_B_new = measure_observable_matrix(rho_B, op, projective=projective_mode, gamma_value)
            
            dict_A[lbl][k] = val_A
            dict_B[lbl][k] = val_B
            
            # Solo sumamos a la distancia si es 1-body o 2-body estándar (para no distorsionar tu métrica original)
            # o podrías sumar todo, depende de tu criterio. Aquí mantengo lo estándar para separación.
            if length(lbl) == N # Simple check
                d_sq += (val_A - val_B)^2
            end
            
            if projective_mode; rho_A, rho_B = rho_A_new, rho_B_new; end
        end
        separation_dist[k] = sqrt(d_sq)
        if k%50==0; print("\rStep $k / $steps"); end
    end

    # --- POST-PROCESADO Y STM ---
    println("\n🧠 Calculando Capacidad de Memoria...")
    capacities, total_stm = calculate_stm_capacity(dict_B, inputs, MAX_DELAY, WASHOUT)
    plot_memory_capacity(capacities, MAX_DELAY, total_stm, Experiment_path)

    # --- VISUALIZACIÓN ---
    println("📊 Generando Gráficas...")
    
    # 1. Plots Estándar (Tus boxplots y validaciones normales)
    for basis in ["X", "Y", "Z"]
        basis_dir = joinpath(Experiment_path, "Basis_$basis"); mkpath(basis_dir)
        sub_A = filter(p -> contains(p.first, basis), dict_A)
        sub_B = filter(p -> contains(p.first, basis), dict_B)
        
        plot_separability_boxplots(sub_A, inputs, N, basis_dir)
        plot_quick_validation_per_qubit(separation_dist, sub_B, inputs, N, basis_dir)
        try; plot_and_save_validation_full(sub_A, sub_B, separation_dist, N, steps, basis_dir); catch; end
    end

    # ==========================================================
    # 2. COMPARATIVA EXHAUSTIVA CON NATHAN (TODAS LAS OPCIONES)
    # ==========================================================
    println("\n⚔️  Iniciando Comparativa Masiva con Nathan...")
    comp_dir = joinpath(Experiment_path, "Comparativa_Nathan_FULL")
    if !isdir(comp_dir); mkpath(comp_dir); end
    
    # Identificamos todas las claves comunes
    common_keys = intersect(keys(nathan_results), keys(dict_B))
    
    if isempty(common_keys)
        println("❌ No hay claves comunes. Revisa los nombres.")
    else
        println("✅ Encontradas $(length(common_keys)) coincidencias. Generando plots...")
        
        # Array para guardar errores
        errors = Float64[]
        labels_err = String[]

        for key in common_keys
            # Plot individual
            p_comp = plot(title="Comparativa: $key", xlabel="Steps", ylabel="<O>", legend=:topright)
            plot!(p_comp, nathan_results[key], label="Nathan", lw=2.5, color=:black, alpha=0.5)
            plot!(p_comp, dict_B[key], label="Mi Simulación", lw=1.5, color=:red, linestyle=:dash)
            
            # Ajuste de escala visual (Importante para evitar el zoom al ruido)
            ylims!(p_comp, -1.1, 1.1)
            
            savefig(p_comp, joinpath(comp_dir, "Comp_$key.png"))
            
            # Error
            mse = mean((nathan_results[key] .- dict_B[key]).^2)
            push!(errors, mse)
            push!(labels_err, key)
            
            # Print discreto si el error es grande
            if mse > 1e-3
                println("   ⚠️ Discrepancia en $key: MSE = $mse")
            end
        end
        
        # Plot Resumen de Errores (Barra de errores)
        if !isempty(errors)
            p_err = bar(labels_err, errors, title="MSE por Observable", xrotation=45, legend=false, ylabel="Mean Squared Error")
            savefig(p_err, joinpath(comp_dir, "_Resumen_Errores_MSE.png"))
            println("   📄 Resumen de errores guardado en _Resumen_Errores_MSE.png")
            println("   🌟 Error promedio global: $(mean(errors))")
        end
    end

    save_qrc_results_jld2(N, steps, T_evol, h_val, inputs, dict_B, Experiment_name)
    println("\n✅ Finalizado. Total STM: $total_stm")
end

run_schrodinger_esp_task()