using LinearAlgebra
using Random
using Plots

# ==============================================================================
# 1. CARGA DE MÓDULOS
# ==============================================================================
# Asegúrate de que las rutas son correctas en tu proyecto
include("src/operator_terms/pauli.jl")      
include("src/operator_terms/hamiltonian.jl") # Modificado
include("src/utils/pauli_algebra.jl")       
include("src/utils/dynamics.jl")
include("src/utils/injection.jl")
include("src/utils/initial_state.jl")           # Nuevo (Estado Inicial)
# Importamos el ploteador de correlaciones
include("src/visualization/plot_spin_dynamics.jl")   

# ==============================================================================
# 2. DEFINICIÓN DEL HAMILTONIANO (Correcto: Eq 32 del Paper)
# ==============================================================================
function build_paper_hamiltonian_corrected(n_qubits::Int, J_scale::Float64, h::Float64; seed::Int=1234)
    Random.seed!(seed)
    H = Operator()
    
    # Interacción: + sum J Z Z (J aleatorio [-0.5, 0.5])
    for i in 0:(n_qubits-2)
        J_val = J_scale * (rand() - 0.5) 
        p_zz = PauliString(0, (1 << i) | (1 << (i+1)))
        H[p_zz] = get(H, p_zz, 0.0im) + J_val
    end
    
    # Campo Transversal: - h sum X (El signo menos es vital)
    for i in 0:(n_qubits-1)
        p_x = PauliString(1 << i, 0)
        H[p_x] = get(H, p_x, 0.0im) - h 
    end
    
    return H
end

# ==============================================================================
# 3. TAREA DE MEMORIA (Short Term Memory)
# ==============================================================================
function main_memory_task()
    println("=== SHORT TERM MEMORY TEST (Random Inputs) ===")
    
    # --- A. CONFIGURACIÓN ---
    n_qubits = 6   # Tamaño del sistema
    n_steps = 100  # Duración del test
    
    # Física (Parámetros ajustados para memoria)
    Js = 1.0
    h = 1.5        # Campo fuerte para propagar la señal
    dt = 0.05      # Paso de tiempo entre inputs
    g = 0.1        # Ruido (Dephasing)
    
    sub_steps = 50 # Pasos de integración RK4 por cada paso de tiempo
    
    # --- B. INPUTS ALEATORIOS ---
    # Generamos una secuencia binaria aleatoria {0, 1}
    # [cite_start]Esto cumple con la nota: "Generate signal list {sk} of 0s, 1s" [cite: 93]
    inputs_sk = rand([0.0, 1.0], n_steps) 
    
    # --- C. CONSTRUCCIÓN DEL SISTEMA ---
    println("Construyendo Hamiltoniano y Estado Inicial...")
    H_ising = build_paper_hamiltonian_corrected(n_qubits, Js, h; seed=1234)
    
    # Schrödinger: dρ/dt = -i[H, ρ]
    H_evol = Operator()
    for (p, c) in H_ising; H_evol[p] = -c; end
    
    # Estado inicial sparse (X1) o All Zeros
    rho = initial_state_x_first(n_qubits)
    
    # Matriz para guardar la magnetización Z de cada qubit
    results_Z = zeros(Float64, n_steps, n_qubits)
    
    println("Simulando secuencia de memoria (Inject -> Dephase -> Evolve)...")
    
    for k in 1:n_steps
        s_k = inputs_sk[k]
        
        # --- PASO 1: CODIFICACIÓN E INYECCIÓN (Inject) ---
        # [cite_start]Implementa la fórmula de la nota[cite: 98, 101]:
        # |psi> = sqrt(1-s)|0> + sqrt(s)|1>
        # rz = 1 - 2s (Si s=0 -> rz=1; Si s=1 -> rz=-1)
        # rx = 2*sqrt(s*(1-s)) (Si s es binario puro, rx será 0)
        
        rz_val = 1.0 - 2.0 * s_k
        rx_val = 2.0 * sqrt(s_k * (1.0 - s_k))
        
        rho = inject_state(rho, 0, rz_val, rx=rx_val)
        
        # --- PASO 2: DEPHASING (Measure/Backaction) ---
        # [cite_start]Orden corregido según notas manuscritas (Pág 2 del PDF) [cite: 2]
        rho = apply_global_dephasing(rho, g)
        
        # --- PASO 3: EVOLUCIÓN (Evolve) ---
        dt_small = dt / sub_steps
        for _ in 1:sub_steps
            rho = step_rk4(rho, H_evol, dt_small)
            # Truncamiento suave para mantener eficiencia
            truncate_operator!(rho, 2000)
        end
        
        # --- PASO 4: MEDIDA (Readout) ---
        # Leemos Z en todos los qubits para ver dónde está la información
        for q in 0:(n_qubits-1)
            p_meas = PauliString(0, 1 << q) # Medir Z
            val = get(rho, p_meas, 0.0im)
            results_Z[k, q+1] = real(val)
        end
    end
    
    println("Terminado.")
    
    # --- VISUALIZACIÓN ---
    println("Generando gráfico de memoria...")
    
    # Graficamos Input vs Output
    # Qubit 1 (Input) debería seguir a s_k
    # Qubit 2, 3... deberían seguirlo con retraso (Memoria)
    p = plot(results_Z, 
        label = reshape(["Q$i" for i in 0:n_qubits-1], 1, n_qubits),
        title = "Short Term Memory Test (Random Inputs)",
        xlabel = "Time Step",
        ylabel = "Magnetization <Z>",
        lw = 1.5,
        layout = (n_qubits, 1), # Un subplot por qubit para verlo claro
        size = (800, 1000),
        legend = :outertopright
    )
    
    # Guardar
    savefig(p, "memory_test_results.png")
    println("Gráfico guardado en 'memory_test_results.png'")
end

main_memory_task()