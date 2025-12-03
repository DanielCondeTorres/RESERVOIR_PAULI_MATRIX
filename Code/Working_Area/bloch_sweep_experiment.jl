using LinearAlgebra
using Random

# --- CARGAR MÓDULOS ---
include("src/operator_terms/pauli.jl")
include("src/operator_terms/hamiltonian.jl")
include("src/utils/pauli_algebra.jl")
include("src/utils/dynamics.jl")
include("src/utils/injection.jl")
# --- NUEVOS MÓDULOS NECESARIOS ---
include("src/utils/initial_state.jl")           # Nuevo (Estado Inicial)

include("src/visualization/expectation_xj_vs_step.jl") 

function main()
    println("=== REPLICANDO FIGURAS 4 y 5 DEL PAPER (Corregido) ===")
    
    # 1. Configuración del Sistema (N=6)
    n_qubits = 20
    k_max = 50
    
    # --- PARÁMETROS DE TRUNCACIÓN ---
    MAX_PAULI_STRINGS = 1024  
    
    # --- FÍSICA DEL PAPER ---
    Js = 1.0
    h_const = 10.0 * Js   
    dt = 10.0 / Js        
    g_strength = 0.3      
    
    println("Params: N=$n_qubits, MaxStrings=$MAX_PAULI_STRINGS, h=$h_const, dt=$dt, g=$g_strength")
    
    # 2. Hamiltoniano (Usamos la nueva función con distribución uniforme) 
    # Nota: build_paper_hamiltonian ya incluye Random.seed! internamente
    H_ising = build_paper_hamiltonian(n_qubits, Js, h_const; seed=1234)
    
    # Preparamos para ecuación de Schrödinger: dρ/dt = -i[H, ρ] [cite: 37]
    H_evol = Operator()
    for (p, c) in H_ising; H_evol[p] = -c; end
    
    # 3. Estado Inicial |000...0> 
    # Antes tenías rho = Operator(), que está vacío. Esto es incorrecto para rho.
    println("Inicializando estado |000...0>...")
    rho = initial_state_all_zeros(n_qubits)
    
    results_X = zeros(Float64, k_max + 1, n_qubits)
    steps_vec = 0:k_max
    
    println("Iniciando simulación...")
    
    for k in 0:k_max
        # A) Vector de Bloch (Sweep adiabático del Paper)
        ratio = k / k_max
        current_rx = sqrt(ratio)       
        current_rz = sqrt(1.0 - ratio) 
        
        # B) INYECCIÓN (Inject)
        # Reseteamos Q0 y escribimos el nuevo estado
        rho = inject_state(rho, 0, current_rz, rx=current_rx)
        
        # C) EVOLUCIÓN TEMPORAL (Evolve) [cite: 6]
        # Primero la dinámica unitaria (Hamiltoniano)
        sub_steps = 2500 
        dt_small = dt / sub_steps
        
        for _ in 1:sub_steps
            rho = step_rk4(rho, H_evol, dt_small)
            truncate_operator!(rho, MAX_PAULI_STRINGS)
        end
        
        # D) RUIDO (Dephasing)
        # Aplicamos el ruido después de la evolución (o como parte de ella)
        # Usamos la nueva función exponencial definida en quantum_channels.jl
        rho = apply_global_dephasing(rho, g_strength)
        
        # E) MEDICIÓN <X> (Measure) [cite: 7]
        # Medimos al final del ciclo, justo antes del siguiente paso
        for q in 0:(n_qubits-1)
            p_meas = PauliString(1 << q, 0) # Medir X
            val = get(rho, p_meas, 0.0im)
            results_X[k+1, q+1] = real(val)
        end
        
        if k % 5 == 0; print("."); end
    end
    println("\nTerminado.")
    
    plot_expectation_evolution(steps_vec, results_X, n_qubits)
end

main()