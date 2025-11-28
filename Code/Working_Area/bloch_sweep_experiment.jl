using LinearAlgebra
using Random

# --- CARGAR MÓDULOS ---
include("src/operator_terms/hamiltonian.jl")
include("src/utils/injection.jl")
include("src/utils/dynamics.jl")
include("src/visualization/expectation_xj_vs_step.jl") # Asegúrate que este archivo exista

function main()
    println("=== REPLICANDO FIGURAS 4 y 5 DEL PAPER ===")
    
    # 1. Configuración del Sistema (N=6)
    n_qubits = 6
    k_max = 50
    
    # --- PARÁMETROS DE TRUNCACIÓN (Ajustar según qué figura quieras) ---
    # Figura 3 (Perfecta): 5000
    # Figura 4 (Half):     2048
    # Figura 5 (Quarter):  1024
    MAX_PAULI_STRINGS = 1024  
    
    # --- FÍSICA DEL PAPER (Sección 3.1) ---
    Js = 1.0
    h_const = 10.0 * Js   # Campo FUERTE (esto causaba tus líneas verticales)
    dt = 10.0 / Js        # Tiempo LARGO
    g_strength = 0.3      # Ruido Dephasing
    
    println("Params: N=$n_qubits, MaxStrings=$MAX_PAULI_STRINGS, h=$h_const, dt=$dt, g=$g_strength")
    
    # Semilla fija para reproducir el desorden
    Random.seed!(1234)

    # 2. Hamiltoniano (J random, h constante)
    H_ising = build_paper_hamiltonian(n_qubits, Js, h_const)
    H_evol = Operator()
    for (p, c) in H_ising; H_evol[p] = -c; end
    
    rho = Operator()
    results_X = zeros(Float64, k_max + 1, n_qubits)
    steps_vec = 0:k_max
    
    println("Iniciando simulación (puede tardar por los sub-pasos)...")
    
    for k in 0:k_max
        # A) Vector de Bloch
        # Para obtener la curva QUE SUBE (Figuras), invertimos el orden del paper
        # rx va de 0 a 1.
        ratio = k / k_max
        rx = sqrt(ratio)
        rz = sqrt(1.0 - ratio)
        
        # B) Inyección
        rho = inject_general_bloch_state(rho, 0, rx, 0.0, rz)
        
        # C) Dephasing (Ruido)
        rho = apply_global_dephasing(rho, g_strength)
        
        # D) Medición <X>
        for q in 0:(n_qubits-1)
            p_meas = PauliString(1 << q, 0)
            val = get(rho, p_meas, 0.0im)
            results_X[k+1, q+1] = real(val)
        end
        
        # E) EVOLUCIÓN ESTABILIZADA
        # Como h=10 y dt=10, el cambio total es enorme (~100).
        # RK4 necesita pasos pequeños (dt ~ 0.005) para no explotar.
        # dt / sub_steps < 0.005 => 10 / sub_steps < 0.005 => sub > 2000
        sub_steps = 2500 
        dt_small = dt / sub_steps
        
        for _ in 1:sub_steps
            rho = step_rk4(rho, H_evol, dt_small)
            # Truncamos frecuentemente para mantener la memoria baja
            truncate_operator!(rho, MAX_PAULI_STRINGS)
        end
        
        if k % 5 == 0; print("."); end
    end
    println("\nTerminado.")
    
    plot_expectation_evolution(steps_vec, results_X, n_qubits)
end

main()