using LinearAlgebra

# --- IMPORTACIÓN DE MÓDULOS ---
include("src/operator_terms/pauli.jl")      
include("src/operator_terms/hamiltonian.jl")
include("src/utils/pauli_algebra.jl")       
include("src/utils/dynamics.jl")
include("src/utils/injection.jl")           

# IMPORTAMOS LOS DOS PLOTEADORES
include("src/visualization/plot_spin_dynamics.jl")      
include("src/visualization/expectation_xj_vs_step.jl")  

function run_correlations_experiment()
    # --- CONFIGURACIÓN ---
    n_qubits = 6
    J = 1.0     
    h = 1.5      # Campo fuerte necesario para propagar la correlación
    dt = 0.05   
    gamma_dephasing = 0.0001 
    
    steps_per_input = 200 
    n_inputs = 100
    
    println("=== QRC Protocol: Correlaciones con el Input (Z0Zk) ===")
    
    # --- HAMILTONIANO ---
    H_ising = build_ising_hamiltonian(n_qubits, J, h)
    H_schrodinger = Operator()
    for (p, c) in H_ising; H_schrodinger[p] = -c; end
    
    # --- SIMULACIÓN ---
    inputs = [isodd(i) ? 0.9 : -0.9 for i in 1:n_inputs]
    rho = Operator()
    rho = inject_state(rho, 0, inputs[1]) 
    
    # --- DEFINICIÓN DE COLUMNAS ---
    # Queremos: [Z0...Z5] seguido de [Z0Z1...Z0Z5]
    # N sitios + (N-1) correlaciones
    n_correlations = n_qubits - 1
    total_cols = n_qubits + n_correlations
    
    times = Float64[]
    results_ordered = zeros(Float64, n_inputs, total_cols)
    
    # 1. Generar Etiquetas en el orden solicitado
    labels_list = String[]
    
    # Bloque A: Sitios individuales
    for i in 0:(n_qubits-1)
        push!(labels_list, "Z$i")
    end
    
    # Bloque B: Correlaciones con Z0
    for i in 1:(n_qubits-1)
        push!(labels_list, "Z0Z$i")
    end
    
    labels_matrix = reshape(labels_list, 1, total_cols)
    println("Observables a medir: $labels_list")
    
    println("Simulando...")
    
    for k in 1:n_inputs
        # A) INYECCIÓN
        if k > 1; rho = inject_state(rho, 0, inputs[k]); end
        
        # B) EVOLUCIÓN
        for s in 1:steps_per_input
            rho = step_rk4(rho, H_schrodinger, dt)
            rho = apply_global_dephasing(rho, gamma_dephasing)
        end
        
        # C) MEDICIÓN (EN EL ORDEN SOLICITADO)
        push!(times, k)
        
        col_idx = 1
        
        # --- PARTE 1: MEDIR SITIOS (Z0, Z1, Z2...) ---
        for q in 0:(n_qubits-1)
            # Pauli Z en posición q
            p_site = PauliString(0, 1 << q)
            val = real(get(rho, p_site, 0.0im))
            results_ordered[k, col_idx] = val
            col_idx += 1
        end
        
        # --- PARTE 2: MEDIR CORRELACIONES CON INPUT (Z0Z1, Z0Z2...) ---
        # El input siempre es Z0 (bit 0)
        bit_0 = 1 << 0
        
        for q in 1:(n_qubits-1)
            bit_q = 1 << q
            # Máscara combinada: Bit 0 Y Bit q encendidos
            mask_corr = bit_0 | bit_q
            
            p_corr = PauliString(0, mask_corr)
            val = real(get(rho, p_corr, 0.0im))
            results_ordered[k, col_idx] = val
            col_idx += 1
        end
    end
    
    println("Simulación terminada.")

    # --- GENERACIÓN DE GRÁFICOS ---
    
    # 1. GRÁFICO COMPLETO (Sitios + Correlaciones Z0Zk)
    # Llama a plot_spin_dynamics.jl
    plot_mixed_dynamics(times, results_ordered, labels_matrix)
    
    # 2. GRÁFICO SOLO SITIOS
    # Extraemos solo las primeras N columnas (Z0...Z5)
    results_sites_only = results_ordered[:, 1:n_qubits]
    
    # Llama a expectation_xj_vs_step.jl
    plot_expectation_evolution(times, results_sites_only, n_qubits)
    
    println("✅ Gráficos generados.")
end

run_correlations_experiment()