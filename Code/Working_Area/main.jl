using LinearAlgebra
using Statistics
using Plots
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
include("src/utils/measurements.jl")

# Feature Extraction
include("src/capacity_training/feature_extraction.jl") 

# Training & Plotting
include("src/capacity_training/qrc_training.jl")    
include("src/visualization/plot_stm_capacity.jl") 

# 🔴 NUEVO: Incluimos tu función de graficar expectación
include("src/visualization/expectation_xj_vs_step.jl") 

# ==============================================================================
# 2. PARÁMETROS DE ENTRADA
# ==============================================================================
N_QUBITS    = 6             
NUM_STEPS   = 1000 #       
DT          = 0.01          
N_SUBSTEPS  = 100           
J_SCALE     = 1.0           
H_FIELD     = 10.0          
G_DEPHASING = 0.05          
PAULI_BASIS = "Z"           
NOISE_SHOTS = 1.5e6         
MAX_DELAY   = 10            
WASHOUT     = 200           

# ==============================================================================
# 3. SCRIPT PRINCIPAL
# ==============================================================================
function run_simulation_from_scratch()
    println("🚀 INICIANDO SIMULACIÓN COMPLETA (STM + Expectation Plot)...")

    # A. Generar Datos Aleatorios
    s_vec = rand([0.0, 1.0], NUM_STEPS)

    # B. Preparar Física
    # Generar desorden J en [-0.5, 0.5] * J_SCALE
    Jvec = (rand(N_QUBITS) .- 0.5) .* J_SCALE
    H_evol = hamiltonian_nathan_XX(N_QUBITS, Jvec, H_FIELD)
    
    rho = initial_state_all_zeros(N_QUBITS)
    dt_small = DT / N_SUBSTEPS
    scale_factor = 2.0^N_QUBITS # 🔴 IMPORTANTE: Factor 2^N

    # C. Preparar Reservorio
    basis = build_reservoir_basis(N_QUBITS)
    n_features = length(basis)
    
    # Matriz para STM (Features)
    X_reservoir = zeros(Float64, NUM_STEPS, n_features)

    # 🔴 NUEVO: Matriz para guardar la evolución de <X_j>
    # Dimensiones: [Pasos x Qubits] -> [1000 x 6]
    X_magnetization_data = zeros(Float64, NUM_STEPS, N_QUBITS)

    # D. BUCLE DE SIMULACIÓN
    println("⏳ Ejecutando evolución ($NUM_STEPS pasos)...")
    
    for k in 1:NUM_STEPS
        # 1. Inyección
        s_k = s_vec[k]
        rz = 1.0 - 2.0 * s_k
        rx = 2.0 * sqrt(s_k * (1.0 - s_k))
        rho = inject_state_EraseWrite(rho, 0, rz, rx=rx)

        # 2. Evolución
        for _ in 1:N_SUBSTEPS
            rho = step_rk4(rho, H_evol, dt_small)
            truncate_operator!(rho, 2000)
        end

        # 3. Medición (Features y Expectation)
        feats_vals = extract_all_features(rho, basis)
        
        # --- Guardar Features para STM (con ruido y escala) ---
        for i in 1:n_features
            val_fisico = feats_vals[i] * scale_factor
            X_reservoir[k, i] = apply_shot_noise(val_fisico, NOISE_SHOTS)
        end

        # 🔴 NUEVO: Guardar datos específicos de <X_j> para la gráfica
        # Sabemos que en 'basis', primero van los Z (N), luego los X (N).
        # Por tanto, los X están en los índices [N+1] hasta [2N].
        offset_X = N_QUBITS 
        for qubit_idx in 1:N_QUBITS
            # Extraemos el valor de la base correspondiente a X_j
            raw_val = feats_vals[offset_X + qubit_idx] 
            
            # Escalamos a valor físico real
            phys_val = raw_val * scale_factor
            
            # Guardamos en la matriz de magnetización
            X_magnetization_data[k, qubit_idx] = phys_val
        end

        # 4. Backaction (g * dt)
        g_effective = G_DEPHASING * DT
        rho = apply_global_dephasing(rho, g_effective, PAULI_BASIS)
        
        if k % 100 == 0; print("\r   -> Progreso: $(round(k/NUM_STEPS*100, digits=1))%"); end
    end
    println("\n✅ Simulación completada.")

    # --------------------------------------------------------------------------
    # E. GRAFICAR EVOLUCIÓN DE EXPECTACIÓN (LO QUE PEDISTE)
    # --------------------------------------------------------------------------
    # Llamamos a tu función pasando: Eje X, Matriz de Datos X_j, Num Qubits
    plot_expectation_evolution(1:NUM_STEPS, X_magnetization_data, N_QUBITS)
    
    # --------------------------------------------------------------------------
    # F. CALCULAR Y GRAFICAR CAPACIDAD STM
    # --------------------------------------------------------------------------
    println("📉 Calculando STM Capacity...")
    capacities = Float64[]

    for tau in 1:MAX_DELAY
        # (Lógica estándar de capacidad...)
        y_target = s_vec[1 : end-tau]
        X_feats  = X_reservoir[1+tau : end, :]
        len_data = length(y_target)
        split_idx = floor(Int, len_data * 0.8)
        
        if split_idx <= WASHOUT; push!(capacities, 0.0); continue; end

        X_train = X_feats[1:split_idx, :]
        y_train = y_target[1:split_idx]
        X_test  = X_feats[split_idx+1:end, :]
        y_test  = y_target[split_idx+1:end]
        
        weights, _ = train_reservoir(X_train, y_train, WASHOUT)
        rows_test = size(X_test, 1)
        X_test_bias = hcat(X_test, ones(rows_test))
        y_pred_test = X_test_bias * weights
        
        C = calculate_capacity(y_test, y_pred_test)
        push!(capacities, C)
        println("   Tau $tau -> C = $(round(C, digits=4))")
    end
    
    total_stm = sum(capacities)
    println("🏆 Capacidad Total STM = $(round(total_stm, digits=4))")
    plot_stm_capacity(capacities, MAX_DELAY, N_QUBITS)
end

run_simulation_from_scratch()