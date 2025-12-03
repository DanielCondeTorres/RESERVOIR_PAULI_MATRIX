using LinearAlgebra
using Random

# --- CARGAR MÓDULOS ---
include("src/operator_terms/pauli.jl")
include("src/operator_terms/hamiltonian.jl")
include("src/utils/pauli_algebra.jl")
include("src/utils/dynamics.jl")
include("src/utils/injection.jl")
include("src/utils/initial_state.jl")           
include("src/visualization/expectation_xj_vs_step.jl") 

function main()
    println("=== REPLICANDO FIGURA 2 DEL PAPER (N=20, Rápido) ===")
    
    # 1. Configuración
    n_qubits = 6
    k_max = 50
    
    # TRUNCAMIENTO AGRESIVO (Para Fig 2) [cite: 266]
    # Usamos 100 como pusiste en tu código anterior
    HARD_LIMIT_TERMS = 800
    
    # Física
    Js = 1.0
    h_const = 10.0 * Js   
    dt = 10.0 / Js        
    g_strength = 0.3      
    
    println("Params: N=$n_qubits, Truncation=$HARD_LIMIT_TERMS")
    
    # 2. Hamiltoniano
    H_ising = build_paper_hamiltonian(n_qubits, Js, h_const; seed=1234)
    H_evol = Operator()
    for (p, c) in H_ising; H_evol[p] = -c; end
    
    # 3. ESTADO INICIAL EFICIENTE 
    # En lugar de calcular el millón de términos de |00..0>, 
    # empezamos con X en el primer qubit (como dice el paper para N grande).
    println("Inicializando rho = X_0 (1 término)...")
    rho = initial_state_x_first(n_qubits)
    
    results_X = zeros(Float64, k_max + 1, n_qubits)
    steps_vec = 0:k_max
    
    println("Iniciando simulación...")
    
    # OPTIMIZACIÓN: 400 pasos es suficiente estabilidad para este dt
    sub_steps = 3000 
    dt_small = dt / sub_steps
    
    for k in 0:k_max
        # A) Vector de Bloch (Sweep Z -> X)
        ratio = k / k_max
        current_rx = sqrt(ratio)       
        current_rz = sqrt(1.0 - ratio) 
        
        # B) INYECCIÓN (Inject)
        rho = inject_state(rho, 0, current_rz, rx=current_rx)
        
        # C) DEPHASING (Measure/Backaction) 
        # Orden corregido según notas manuscritas [cite: 334]
        rho = apply_global_dephasing(rho, g_strength)
        
        # D) EVOLUCIÓN (Evolve)
        for _ in 1:sub_steps
            rho = step_rk4(rho, H_evol, dt_small)
            
            # Truncamiento constante para velocidad
            #truncate_operator!(rho, HARD_LIMIT_TERMS)
            truncate_operator!(rho)
        end
        
        # E) MEDICIÓN
        for q in 0:(n_qubits-1)
            p_meas = PauliString(1 << q, 0)
            val = get(rho, p_meas, 0.0im)
            results_X[k+1, q+1] = real(val)
        end
        
        print("\rPaso $k / $k_max listo.")
    end
    println("\nTerminado.")
    
    plot_expectation_evolution(steps_vec, results_X, n_qubits)
end

main()