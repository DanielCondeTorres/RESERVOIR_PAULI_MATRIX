using LinearAlgebra
using Random
using Plots

# --- CARGAR MÓDULOS ---
include("src/operator_terms/pauli.jl")      
include("src/operator_terms/hamiltonian.jl")
include("src/utils/pauli_algebra.jl")       
include("src/utils/dynamics.jl")
include("src/utils/injection.jl")
include("src/utils/initial_state.jl")           
include("src/utils/qrc_training.jl") 
include("src/visualization/plot_stm_capacity.jl") 

# IMPORTANTE: Cargamos el nuevo módulo de mediciones
include("src/utils/measurements.jl") 

function run_paper_replica_stm()
    println("=== REPLICANDO FIGURA 3a (STM Capacity) - Código Limpio ===")
    
    # 1. PARÁMETROS
    n_qubits = 6
    n_steps = 1000
    washout = 200
    Js = 1.0; h = 10.0; dt = 10.0; g = 0.3
    
    sub_steps = 1000
    dt_small = dt / sub_steps
    
    # 2. PREPARACIÓN DE MEDICIONES (Esto es nuevo)
    # Generamos la lista de observables (PauliStrings) UNA sola vez fuera del bucle.
    # Esto es mucho más eficiente y limpio.
    observables_list = build_paper_basis(n_qubits)
    n_features = length(observables_list) # 17 automáticamente
    
    # 3. INICIALIZACIÓN
    inputs_sk = rand(n_steps) 
    H_ising = build_paper_hamiltonian(n_qubits, Js, h; seed=1234)
    H_evol = Operator()
    for (p, c) in H_ising; H_evol[p] = -c; end
    rho = initial_state_x_first(n_qubits)
    
    # Matriz para guardar resultados
    reservoir_states = zeros(Float64, n_steps, n_features)
    
    println("Simulando $n_steps pasos con $n_features features...")
    
    for k in 1:n_steps
        s_k = inputs_sk[k]
        
        # A) INYECCIÓN
        rz_val = 1.0 - 2.0 * s_k
        rx_val = 2.0 * sqrt(s_k * (1.0 - s_k)) 
        rho = inject_state(rho, 0, rz_val, rx=rx_val)
        
        # B) EVOLUCIÓN
        for _ in 1:sub_steps
            rho = step_rk4(rho, H_evol, dt_small)
            truncate_operator!(rho, 2000) 
        end
        
        # C) DEPHASING (Ruido)
        rho = apply_global_dephasing(rho, g)
        
        # D) READOUT (¡AHORA ES UNA SOLA LÍNEA!)
        # Usamos la función extract_features que definimos arriba.
        # projective=false porque para la Figura 3a queremos los valores esperados (Expectation Values)
        # tal como se comportaría un ensemble infinito (o promediado).
        features_vec = extract_features(rho, observables_list)
        
        # Guardamos el vector en la fila k
        reservoir_states[k, :] = features_vec
        
        if k % 50 == 0; print("\rPaso $k..."); end
    end
    println("\nSimulación finalizada.")

    # 4. ENTRENAMIENTO
    println("Calculando capacidades...")
    max_delay = 10
    capacities = zeros(Float64, max_delay)
    
    for tau in 1:max_delay
        target = inputs_sk[1:end-tau]
        feats = reservoir_states[tau+1:end, :]
        
        weights, preds = train_reservoir(feats, target, washout)
        
        real_target = target[washout+1:end]
        C = calculate_capacity(real_target, preds[washout+1:end])
        capacities[tau] = C
        println("Tau $tau: C = $C")
    end
    
    plot_stm_capacity(capacities, max_delay, n_qubits)
end

run_paper_replica_stm()