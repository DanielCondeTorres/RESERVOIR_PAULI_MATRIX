using LinearAlgebra
using Random

# --- CARGAR MÓDULOS ---
include("src/operator_terms/hamiltonian.jl")
include("src/utils/injection.jl") # Asegúrate de que inject_state está aquí
include("src/utils/dynamics.jl")
include("src/utils/pauli_algebra.jl")

include("src/visualization/expectation_xj_vs_step.jl") 

function main()
    println("=== REPLICANDO FIGURAS 4 y 5 DEL PAPER (Actualizado) ===")
    
    # 1. Configuración del Sistema (N=6)
    n_qubits = 6
    k_max = 50
    
    # --- PARÁMETROS DE TRUNCACIÓN ---
    # Figura 5 (Quarter):  1024
    MAX_PAULI_STRINGS = 1024  
    
    # --- FÍSICA DEL PAPER ---
    Js = 1.0
    h_const = 10.0 * Js   
    dt = 10.0 / Js        
    g_strength = 0.3      
    
    println("Params: N=$n_qubits, MaxStrings=$MAX_PAULI_STRINGS, h=$h_const, dt=$dt, g=$g_strength")
    
    Random.seed!(1234)

    # 2. Hamiltoniano
    H_ising = build_paper_hamiltonian(n_qubits, Js, h_const)
    H_evol = Operator()
    for (p, c) in H_ising; H_evol[p] = -c; end
    
    rho = Operator()
    results_X = zeros(Float64, k_max + 1, n_qubits)
    steps_vec = 0:k_max
    
    println("Iniciando simulación...")
    
    for k in 0:k_max
        # A) Vector de Bloch (Lógica del Paper)
        ratio = k / k_max
        current_rx = sqrt(ratio)       # Queremos que X crezca
        current_rz = sqrt(1.0 - ratio) # Queremos que Z decrezca
        
        # B) INYECCIÓN CON LA NUEVA FUNCIÓN
        # Firma: inject_state(rho, qubit_idx, rz; rx=0.0, ry=0.0)
        # Nota: Pasamos rz como argumento posicional y rx como keyword.
        rho = inject_state(rho, 0, current_rz, rx=current_rx)
        
        # C) Dephasing (Ruido)
        rho = apply_global_dephasing(rho, g_strength)
        
        # D) Medición <X>
        for q in 0:(n_qubits-1)
            p_meas = PauliString(1 << q, 0)
            val = get(rho, p_meas, 0.0im)
            results_X[k+1, q+1] = real(val)
        end
        
        # E) EVOLUCIÓN ESTABILIZADA
        sub_steps = 2500 
        dt_small = dt / sub_steps
        
        for _ in 1:sub_steps
            rho = step_rk4(rho, H_evol, dt_small)
            truncate_operator!(rho, MAX_PAULI_STRINGS)
        end
        
        if k % 5 == 0; print("."); end
    end
    println("\nTerminado.")
    
    plot_expectation_evolution(steps_vec, results_X, n_qubits)
end

main()