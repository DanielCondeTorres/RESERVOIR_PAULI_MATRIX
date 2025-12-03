using LinearAlgebra
using Random

# --- IMPORTACIÓN DE MÓDULOS ---
include("src/operator_terms/pauli.jl")      
include("src/operator_terms/hamiltonian.jl") # Modificado
include("src/utils/pauli_algebra.jl")       
include("src/utils/dynamics.jl")
include("src/utils/injection.jl")
include("src/utils/initial_state.jl")           # Nuevo (Estado Inicial)

# Importamos el ploteador de correlaciones
include("src/visualization/plot_spin_dynamics.jl")      

function run_paper_compliant_protocol()
    # --- CONFIGURACIÓN ---
    n_qubits = 6
    
    # Parámetros físicos según notas
    J_s = 1.0           # Escala de referencia para J aleatorio
    h = 1.5             # Campo transversal
    dt = 0.05           # Paso de tiempo de integración
    g_strength = 0.1    # Parámetro 'g' para la nueva fórmula de dephasing
    
    steps_per_input = 50 # Pasos de evolución unitaria entre inyecciones
    n_inputs = 100
    
    println("=== QRC Protocol: Cumpliendo Notas y Nueva Fórmula ===")
    println("Params: J_s=$J_s (random), h=$h, g=$g_strength, dt=$dt")
    
    # --- 1. HAMILTONIANO (Nota 38-41) ---
    # Usamos el generador aleatorio según el paper
    H_ising = build_paper_hamiltonian(n_qubits, J_s, h; seed=42)
    
    # Preparamos para Schrödinger: dρ/dt = -i[H, ρ]
    # (Nota 36-37: verificar signo. Tu derivative usa i[H,O], así que pasamos -H)
    H_schrodinger = Operator()
    for (p, c) in H_ising; H_schrodinger[p] = -c; end
    
    # --- 2. INPUTS (Probabilidades) ---
    # Los inputs s_k deben ser probabilidades [0,1] para la fórmula de Bloch
    # Usamos una onda cuadrada de probabilidades (ej: 0.1 y 0.9)
    inputs_sk = [isodd(i) ? 0.9 : 0.1 for i in 1:n_inputs]
    
    # --- 3. ESTADO INICIAL (Nota 46-47) ---
    # Empezamos en |000...0>
    println("Inicializando sistema en estado |000...0>...")
    rho = initial_state_all_zeros(n_qubits)
    
    # --- PREPARACIÓN DE DATOS PARA PLOT ---
    # Mediremos Correlaciones con el Input: Z0Z1, Z0Z2...
    n_correlations = n_qubits - 1
    times = Float64[]
    results = zeros(Float64, n_inputs, n_correlations)
    
    labels_list = ["Z0Z$i" for i in 1:(n_qubits-1)]
    labels_matrix = reshape(labels_list, 1, n_correlations)
    
    println("Simulando ciclo: Inject -> Evolve -> Measure...")
    
    for k in 1:n_inputs
        s_k = inputs_sk[k]
        
        # --- PASO A: INYECCIÓN (Inject) ---
        # Calcular vector de Bloch basado en s_k (Nota 98-101)
        # |ψ_k> = sqrt(1-s_k)|0> + sqrt(s_k)|1>
        # r_z = 1 - 2*s_k
        # r_x = 2 * sqrt(s_k * (1 - s_k))
        rz_val = 1.0 - 2.0 * s_k
        rx_val = 2.0 * sqrt(s_k * (1.0 - s_k))
        
        # Usamos la función unificada inyectando en Z y X
        rho = inject_state(rho, 0, rz_val; rx=rx_val)
        
        # --- PASO B: EVOLUCIÓN (Evolve) ---
        # 1. Evolución Unitaria (Hamiltoniano)
        for s in 1:steps_per_input
            rho = step_rk4(rho, H_schrodinger, dt)
        end
        
        # 2. Evolución No-Unitaria (Dephasing) - APLICADA UNA VEZ POR INPUT
        # Usamos la nueva función con la fórmula exponencial
        rho = apply_global_dephasing(rho, g_strength)
        
        # --- PASO C: MEDICIÓN (Measure) ---
        push!(times, k)
        
        # Medir correlaciones Z0Zi
        bit_0 = 1 << 0
        col_idx = 1
        for q in 1:(n_qubits-1)
            bit_q = 1 << q
            mask_corr = bit_0 | bit_q
            p_corr = PauliString(0, mask_corr)
            val = real(get(rho, p_corr, 0.0im))
            results[k, col_idx] = val
            col_idx += 1
        end
    end
    
    println("Simulación terminada.")
    println("Últimas correlaciones Z0Zi: $(round.(results[end, :], digits=4))")

    # --- GRAFICAR ---
    # Usamos el ploteador mixto para ver las correlaciones
    plot_mixed_dynamics(times, results, labels_matrix)
end

run_paper_compliant_protocol()