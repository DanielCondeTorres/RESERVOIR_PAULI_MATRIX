using LinearAlgebra
using Statistics
using Plots
using JLD2
using Random  # <--- ESTO ARREGLA TU ERROR ACTUAL
using Printf

# ==============================================================================
# 1. CARGA DE MÓDULOS (Asegúrate que estas rutas existen en tu carpeta src/)
# ==============================================================================
# Ajusta las rutas si es necesario. Asumo que corres desde la carpeta raíz.
include("src/operator_terms/pauli_algebra.jl")      
include("src/utils/dynamics.jl")
include("src/utils/injection_EraseWrite.jl")
include("src/utils/initial_state.jl") 
include("src/utils/quantum_channels.jl") 
include("src/utils/shot_noise.jl") 
include("src/operator_terms/hamiltonian.jl") 
include("src/loaders/data_loader.jl")
include("src/utils/measurements.jl")
# Si tienes estas, inclúyelas. Si no, las funciones abajo suplen la parte de cálculo.
include("src/capacity_training/feature_extraction.jl") 
# include("src/capacity_training/qrc_training.jl") <--- Usaremos la versión Ridge abajo
include("src/visualization/plot_stm_capacity.jl") 



# Variables globales que tu código original esperaba
INPUT_FILE = "6_1_2_all_zeros_12345.jld2" 
N_SUBSTEPS = 50 # Reducido un poco para velocidad, ajusta si quieres más precisión
MAX_DELAY   = 10        
WASHOUT     = 200       
NOISE_SHOTS = 1.5e6     
DEPHASING_PAULI_MATRIX = "Z" 

function run_stm_final()
    println("🚀 Iniciando Simulación de Capacidad STM (FULL READY)...")

    # --- 1. DEFINICIÓN DE PARÁMETROS (Blindada) ---
    params = Dict(
        "N_qubits" => 6, 
        "num_steps" => 3000, 
        "J" => 1.0, 
        "h" => 0.5, 
        "g" => 0.1, 
        "dt" => 0.1
    )

    # Intento de carga de archivo
    if isfile(INPUT_FILE)
        try
            loaded = extract_metadata(INPUT_FILE)
            merge!(params, loaded)
            println("📂 Metadatos leídos de $INPUT_FILE")
        catch
            println("⚠️ Error leyendo archivo. Usando defaults.")
        end
    else
        println("⚠️ Archivo no encontrado. Usando defaults.")
    end

    # --- 2. FORZADO DE RÉGIMEN FÍSICO (CRÍTICO) ---
    N = params["N_qubits"]
    n_steps = 3000   # Pasos suficientes para estadística
    dt_step = 4.0    # TIEMPO LARGO: Permite que los espines interactúen antes de borrar
    #params["J"] = 1.0 
    
    println("🔹 Config: N=$N | Steps=$n_steps | dt=$dt_step | J=1.0 (Régimen de mezcla fuerte)")

    # --- 3. GENERACIÓN DE SEÑAL UNIFICADA ---
    Random.seed!(1234)
    s_vec_unified = rand(Float64, n_steps) # Señal única para inyección y target
    
    # --- 4. PREPARACIÓN FÍSICA ---
    # Generamos J aleatorio pero fuerte
    J_width = params["J"]
    J_vec = (rand(N) .- 0.5) .* J_width 
    J_vec = sign.(J_vec) .* max.(abs.(J_vec), 0.1) # Evita J muy pequeños
    H_evol = hamiltonian_nathan_XX(N, J_vec, params["h"])
    
    rho = initial_state_all_zeros(N)
    
    dt_rk4  = dt_step / N_SUBSTEPS
    g_effective = params["g"] * dt_step 

    basis = build_reservoir_basis(N)
    n_features = length(basis)
    X_reservoir = zeros(Float64, n_steps, n_features)

    println("⏳ Evolucionando sistema...")
    # --- 5. BUCLE TEMPORAL ---
    for k in 1:n_steps
        s_k = s_vec_unified[k]
        
        # A. Inyección (Encoding)
        rz = 1.0 - 2.0 * s_k        
        rx = 2.0 * sqrt(s_k * (1.0 - s_k))
        rho = inject_state_EraseWrite(rho, 0, rz, rx=rx)

        # B. Evolución
        for _ in 1:N_SUBSTEPS
            rho = step_rk4(rho, H_evol, dt_rk4)
        end

        # C. Medición
        raw_coeffs = extract_all_features(rho, basis)
        for i in 1:n_features
            # Clamp en [-1, 1] para que el shot noise funcione
            val = clamp(raw_coeffs[i], -1.0, 1.0)
            X_reservoir[k, i] = apply_shot_noise(val, NOISE_SHOTS)
        end

        # D. Ruido (Dephasing)
        rho = apply_global_dephasing(rho, g_effective, DEPHASING_PAULI_MATRIX)
        
        if k % 500 == 0; print(" $k.."); end
    end
    println("\n✅ Evolución terminada.")

    # --- 6. ENTRENAMIENTO Y CÁLCULO DE C ---
    println("📉 Calculando Capacidad...")
    capacities = Float64[]
    targets_all = s_vec_unified 

    # Revisamos Tau 0 (inmediato) hasta MAX_DELAY
    for tau in 0:MAX_DELAY
        if tau == 0
            y_target = targets_all
            X_feats = X_reservoir
        else
            y_target = targets_all[1 : end-tau]
            X_feats  = X_reservoir[1+tau : end, :]
        end
        
        # Split Train/Test (90% Train para asegurar buen fit)
        len_data = length(y_target)
        split_idx = floor(Int, len_data * 0.9)
        
        if split_idx <= WASHOUT; continue; end

        X_train = X_feats[1:split_idx, :]
        y_train = y_target[1:split_idx]
        X_test  = X_feats[split_idx+1:end, :]
        y_test  = y_target[split_idx+1:end]
        
        # Usamos Ridge Regression
        weights, _ = train_reservoir(X_train, y_train, WASHOUT)
        
        # Predecir
        X_test_bias = hcat(X_test, ones(size(X_test, 1)))
        y_pred = X_test_bias * weights
        
        C = calculate_capacity(y_test, y_pred)
        
        # Solo guardamos para graficar si tau > 0, pero imprimimos todo
        if tau > 0; push!(capacities, C); end
        
        @printf("   Tau %2d -> C = %.4f\n", tau, C)
    end

    total_stm = sum(capacities)
    println("------------------------------------------------")
    println("🏆 Capacidad STM Total (Sum Tau 1-$MAX_DELAY) = $(round(total_stm, digits=4))")
    println("------------------------------------------------")

    # Graficar si es posible
    if isdefined(Main, :plot_stm_capacity)
        plot_stm_capacity(capacities, MAX_DELAY, N)
    elseif isdefined(Main, :Plots)
        # Fallback simple si la función de plot no está cargada
        p = bar(1:length(capacities), capacities, label="STM Capacity", 
                xlabel="Delay (tau)", ylabel="Capacity", title="STM Capacity (Total: $(round(total_stm, digits=2)))")
        display(p)
    end
end

# --- EJECUTAR ---
run_stm_final()