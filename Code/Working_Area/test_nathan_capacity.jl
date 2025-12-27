using LinearAlgebra
using Statistics
using Plots
using JLD2
using Random 

# ==============================================================================
# 1. CARGA DE MÓDULOS
# ==============================================================================
include("src/operator_terms/pauli_algebra.jl")      
include("src/utils/dynamics.jl")
include("src/utils/injection_EraseWrite.jl")
include("src/utils/initial_state.jl") 
include("src/utils/quantum_channels.jl") 
include("src/utils/shot_noise.jl") 
include("src/operator_terms/hamiltonian.jl") 
include("src/loaders/data_loader.jl")
include("src/utils/measurements.jl")

# Feature Extraction & Training
include("src/capacity_training/feature_extraction.jl") 
include("src/capacity_training/qrc_training.jl")             
include("src/visualization/plot_stm_capacity.jl") 

# ==============================================================================
# 2. CONFIGURACIÓN
# ==============================================================================
INPUT_FILE = "6_3_2_all_zeros_12345.jld2" 
N_SUBSTEPS = 100       

# Parámetros STM
MAX_DELAY   = 10        
WASHOUT     = 200       
NOISE_SHOTS = 1.5e6     
DEPHASING_PAULI_MATRIX = "Z" 

# ==============================================================================
# 3. EJECUCIÓN
# ==============================================================================
function run_stm_final()
    println("🚀 Iniciando Simulación de Capacidad STM (Versión Corregida)...")

    # A. Cargar Metadatos
    params = extract_metadata(INPUT_FILE)
    N = params["N_qubits"]
    n_steps = params["num_steps"]
    
    # --- FACTOR DE ESCALA (LA CLAVE QUE FALTABA) ---
    # Convertimos coeficientes matemáticos a valores físicos reales [-1, 1]
    scale_factor = 2.0^N
    println("🔹 Factor de escala aplicado: 2^$N = $scale_factor")

    # B. Preparar Física CON DESORDEN (J Aleatorio)
    J_width = params["J"]
    J_vec = (rand(N) .- 0.5) .* J_width 
    
    # Usamos g del archivo (0.3) pero lo aplicaremos correctamente después
    println("🔹 J_vec generado. Usando g nominal: $(params["g"])")

    H_evol = hamiltonian_nathan_XX(N, J_vec, params["h"])
    rho = initial_state_all_zeros(N)
    
    dt_step = params["dt"] 
    dt_rk4  = dt_step / N_SUBSTEPS

    # C. Preparar Observables
    basis = build_reservoir_basis(N)
    n_features = length(basis)
    println("🔹 Base del reservorio: $n_features observables.")

    X_reservoir = zeros(Float64, n_steps, n_features)

    # D. Bucle Temporal
    println("⏳ Ejecutando evolución...")
    for k in 1:n_steps
        # 1. Injection
        s_k = params["s_vec"][k]
        rz = 1.0 - 2.0 * s_k
        rx = 2.0 * sqrt(s_k * (1.0 - s_k))
        rho = inject_state_EraseWrite(rho, 0, rz, rx=rx)

        # 2. Evolution
        for _ in 1:N_SUBSTEPS
            rho = step_rk4(rho, H_evol, dt_rk4)
            truncate_operator!(rho, 2000)
        end

        # 3. Medición (AQUÍ ESTÁ EL CAMBIO IMPORTANTE)
        raw_coeffs = extract_all_features(rho, basis)
        
        for i in 1:n_features
            # Paso 1: Convertir coeficiente a valor esperado real
            val_fisico = raw_coeffs[i] * scale_factor
            
            # Paso 2: Aplicar ruido al valor físico (que ya tiene tamaño correcto)
            X_reservoir[k, i] = apply_shot_noise(val_fisico, NOISE_SHOTS)
        end

        # 4. Dephasing (Usando g * dt para replicar TASA de decaimiento)
        # Esto hace que g=0.3 se comporte como Nathan quiere (suave) y no como un martillo.
        g_effective = params["g"] * dt_step
        rho = apply_global_dephasing(rho, g_effective, DEPHASING_PAULI_MATRIX)
        
        if k % 100 == 0; print("\r   ⏳ Paso $k / $n_steps ..."); end
    end
    println("\n✅ Datos generados.")

    # Check rápido: Ver si la matriz tiene datos reales o ceros
    avg_val = mean(abs.(X_reservoir))
    println("📊 Valor medio absoluto en el reservorio: $avg_val")
    if avg_val < 1e-4
        println("⚠️ ALERTA: La señal es demasiado débil. Revisa la inyección o el scale_factor.")
    end

    # E. Entrenamiento y Capacidad
    println("📉 Calculando curva de capacidad...")
    capacities = Float64[]
    targets_all = params["s_vec"]

    for tau in 1:MAX_DELAY
        # Alinear Target y Features
        y_target = targets_all[1 : end-tau]
        X_feats  = X_reservoir[1+tau : end, :]
        
        # Split Train/Test
        len_data = length(y_target)
        split_idx = floor(Int, len_data * 0.8)
        
        if split_idx <= WASHOUT
            push!(capacities, 0.0)
            continue
        end

        X_train = X_feats[1:split_idx, :]
        y_train = y_target[1:split_idx]
        X_test  = X_feats[split_idx+1:end, :]
        y_test  = y_target[split_idx+1:end]
        
        # Entrenar
        weights, _ = train_reservoir(X_train, y_train, WASHOUT)
        
        # Predecir
        rows_test = size(X_test, 1)
        X_test_bias = hcat(X_test, ones(rows_test))
        y_pred_test = X_test_bias * weights
        
        # Calcular Capacidad
        C = calculate_capacity(y_test, y_pred_test)
        
        push!(capacities, C)
        println("   Tau $tau -> C = $(round(C, digits=4))")
    end

    total_stm = sum(capacities)
    println("🏆 Capacidad Total STM = $(round(total_stm, digits=4))")

    # F. Graficar
    plot_stm_capacity(capacities, MAX_DELAY, N)
end

run_stm_final()