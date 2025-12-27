using LinearAlgebra
using Statistics
using Plots
using Random # Necesario para generar inputs aleatorios

# ==============================================================================
# 1. CARGA DE MÓDULOS (Tus funciones)
# ==============================================================================
include("src/operator_terms/pauli_algebra.jl")      
include("src/utils/dynamics.jl")
include("src/utils/injection_EraseWrite.jl")
include("src/utils/initial_state.jl") 
include("src/utils/quantum_channels.jl") 
include("src/utils/shot_noise.jl") 
include("src/operator_terms/hamiltonian.jl") 
include("src/utils/measurements.jl")
# Feature Extraction (Base del reservorio)

# --- Extracción de Características (Tu archivo actualizado) ---
include("src/capacity_training/feature_extraction.jl") 

# --- Análisis y Visualización (Tus archivos existentes) ---
include("src/capacity_training/qrc_training.jl")    
include("src/visualization/plot_stm_capacity.jl") 
# --------------------------------------------------------------------------
# A. PARÁMETROS DE ENTRADA (¡MODIFICA ESTO A TU GUSTO!)
# --------------------------------------------------------------------------
# 1. Dimensiones y Tiempo
N_QUBITS    = 6             # Tamaño del sistema
NUM_STEPS   = 1000          # Pasos de tiempo totales
DT          = 0.01          # Tamaño del paso temporal
N_SUBSTEPS  = 100           # Precisión RK4 (substeps por paso)
# 2. Física del Hamiltoniano (Paper: Disordered XX Model)
J_SCALE     = 1.0           # Escala de acoplamiento J
H_FIELD     = 10.0          # Campo magnético transversal h
# 3. Dinámica del Reservorio
G_DEPHASING = 0.05          # Fuerza del ruido (Backaction). 0.05 es ideal.
PAULI_BASIS = "Z"           # Base de dephasing y medición principal
# 4. Realismo Experimental
NOISE_SHOTS = 1.5e6         # Ruido de disparo (infinito = sin ruido)
# 5. Análisis de Memoria
MAX_DELAY   = 10            # Tau máximo a calcular
WASHOUT     = 200           # Pasos iniciales a descartar


# ==============================================================================
# 2. SCRIPT PRINCIPAL (AUTÓNOMO)
# ==============================================================================
function run_simulation_from_scratch()
    println("🚀 INICIANDO SIMULACIÓN QRC (Modo Manual/Generativo)...")
    # --------------------------------------------------------------------------
    # B. GENERACIÓN DE DATOS ALEATORIOS (INPUT SIGNAL)
    # --------------------------------------------------------------------------
    println("🎲 Generando señal de entrada aleatoria...")
    # Generamos un vector de 0s y 1s aleatorios
    s_vec = rand([0.0, 1.0], NUM_STEPS)

    # --------------------------------------------------------------------------
    # C. PREPARACIÓN DE LA FÍSICA
    # --------------------------------------------------------------------------
    println("⚙️ Configurando sistema cuántico...")
    
    # 1. Construir Hamiltoniano XX
    #    Como no cargamos archivo, asumimos que la función genera el desorden internamente
    #    o usa J uniforme. 
    Jvec = (rand(N_QUBITS) .- 0.5) .* J_SCALE
    H_evol = hamiltonian_nathan_XX(N_QUBITS, Jvec, H_FIELD)
    
    # 2. Estado Inicial |00...0>
    rho = initial_state_all_zeros(N_QUBITS)
    
    # 3. Paso de integración pequeño
    dt_small = DT / N_SUBSTEPS

    # 4. Construir Base del Reservorio (Neuronas)
    basis = build_reservoir_basis(N_QUBITS)
    n_features = length(basis)
    println("   -> Reservorio con $N_QUBITS qubits y $n_features observables.")

    # Matriz para guardar la historia [Tiempo x Neuronas]
    X_reservoir = zeros(Float64, NUM_STEPS, n_features)

    # --------------------------------------------------------------------------
    # D. BUCLE DE SIMULACIÓN
    # --------------------------------------------------------------------------
    println("⏳ Ejecutando evolución temporal ($NUM_STEPS pasos)...")
    
    for k in 1:NUM_STEPS
        # 1. INYECCIÓN (Input en Qubit 0)
        s_k = s_vec[k]
        rz = 1.0 - 2.0 * s_k
        rx = 2.0 * sqrt(s_k * (1.0 - s_k))
        rho = inject_state_EraseWrite(rho, 0, rz, rx=rx)

        # 2. EVOLUCIÓN (RK4)
        for _ in 1:N_SUBSTEPS
            rho = step_rk4(rho, H_evol, dt_small)
            truncate_operator!(rho, 2000) # Limpieza numérica
        end

        # 3. MEDICIÓN (Extraer Features)
        feats_vals = extract_all_features(rho, basis)
        
        # Añadir ruido de disparo y guardar
        for i in 1:n_features
            X_reservoir[k, i] = apply_shot_noise(feats_vals[i], NOISE_SHOTS)
        end

        # 4. BACKACTION (Dephasing Global)
        rho = apply_global_dephasing(rho, G_DEPHASING, PAULI_BASIS)
        
        # Progreso visual
        if k % 100 == 0
            print("\r   -> Progreso: $(round(k/NUM_STEPS*100, digits=1))%")
        end
    end
    println("\n✅ Simulación completada.")

    # --------------------------------------------------------------------------
    # E. ENTRENAMIENTO Y CÁLCULO DE CAPACIDAD
    # --------------------------------------------------------------------------
    println("📉 Entrenando Pesos y Calculando STM...")
    capacities = Float64[]

    for tau in 1:MAX_DELAY
        # 1. Definir Target (Lo que pasó hace tau pasos)
        y_target = s_vec[1 : end-tau]
        
        # 2. Definir Features (El estado actual del reservorio)
        X_feats  = X_reservoir[1+tau : end, :]
        
        # 3. Separar Train / Test (80% / 20%)
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
        
        # 4. Entrenar (Regresión Lineal - Ridge/Least Squares)
        weights, _ = train_reservoir(X_train, y_train, WASHOUT)
        
        # 5. Predecir (Test)
        rows_test = size(X_test, 1)
        X_test_bias = hcat(X_test, ones(rows_test)) # Bias manual
        y_pred_test = X_test_bias * weights
        
        # 6. Calcular Capacidad
        C = calculate_capacity(y_test, y_pred_test)
        
        push!(capacities, C)
        println("   Tau $tau -> C = $(round(C, digits=4))")
    end

    total_stm = sum(capacities)
    println("🏆 Capacidad Total STM = $(round(total_stm, digits=4))")

    # --------------------------------------------------------------------------
    # F. GUARDAR RESULTADOS
    # --------------------------------------------------------------------------
    plot_stm_capacity(capacities, MAX_DELAY, N_QUBITS)
    println("💾 Gráfica guardada en la carpeta Outputs.")
end

# Ejecutar script
run_simulation_from_scratch()