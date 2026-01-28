using LinearAlgebra
using Plots
using Random
using JLD2
using FileIO

# ==============================================================================
# 1. CONFIGURACIÓN DE RUTAS Y CARGA
# ==============================================================================
const SCRIPT_DIR = @__DIR__
function include_rel(path...)
    include(joinpath(SCRIPT_DIR, path...))
end

# Cargamos tus archivos
include_rel("src/operator_terms/pauli_algebra.jl")      
include_rel("src/operator_terms/hamiltonian.jl")        
include_rel("src/utils/initial_state.jl") # <--- Tu nueva función está aquí
include_rel("src/utils/measurements.jl")
include_rel("src/utils/dynamics.jl")
include_rel("src/visualization/expectation_xj_vs_step.jl") 

# ==============================================================================
# 2. PARÁMETROS
# ==============================================================================
N = 6
steps = 100
T_evol = 10.0
h_val = 1.0

##### FUNCION AUXILIAR
function dense_to_operator(rho_mat::Matrix{ComplexF64}, n_qubits::Int)
    rho_op = Operator()
    dim = 2^n_qubits
    scale = 1.0 / dim
    # Para medir Z en cada sitio, solo necesitamos los términos de la base Z
    # Pero para ser exactos, convertimos los términos que nos interesan.
    for site in 0:(n_qubits-1)
        z_mask = 1 << site
        p = PauliString(0, z_mask)
        
        # Construimos la matriz del Pauli P para proyectar
        P_mat = [1.0;;]
        for i in 0:(n_qubits-1)
            gate = (i == site) ? [1 0; 0 -1] : [1 0; 0 1]
            P_mat = kron(gate, P_mat)
        end
        
        # Coeficiente en la expansión de Pauli: c_k = Tr(rho * P) / 2^N
        rho_op[p] = tr(rho_mat * P_mat) * scale
    end
    return rho_op
end
# ==============================================================================
# 3. SIMULACIÓN EXACTA
# ==============================================================================
function run_nathan_task()
    println("🚀 Iniciando Simulación Exacta con Matriz de Densidad...")

    # A. Construir Hamiltoniano y Operador de Evolución Exacto
    H_op = build_nathan_all_to_all_XX(N, h_val)
    H_dense = operator_to_dense_matrix(H_op, N)
    U_evol = exp(-im * H_dense * T_evol)
    U_evol_adj = adjoint(U_evol) # U† para la evolución de rho

    # B. INITIAL STATE
    # Operator to density matrix |00...0><00...0|
    println("Create initial state |00...0> ...")
    rho_op = initial_state_all_zeros(N)
    rho = operator_to_dense_matrix(rho_op, N)
   inputs = rand(0:1, steps)
    history = zeros(Float64, steps, N)

    # Pre-calculamos los operadores de medida Z para ganar velocidad
    println("⚙️ Pre-calculando operadores de medida...")
    Z_ops = []
    for site in 0:(N-1)
        Z_site = [1.0;;]
        for i in 0:(N-1)
            gate = (i == site) ? [1 0; 0 -1] : [1 0; 0 1]
            Z_site = kron(gate, Z_site)
        end
        push!(Z_ops, Z_site)
    end

    # D. Bucle de Reservoir
    println("🔄 Ejecutando pasos temporales...")
    for k in 1:steps
        s_k = inputs[k]

        # 1. INYECT STATE
        # Creamos el operador unitario de inyección U_inj
        theta = 2.0 * asin(sqrt(s_k))
        Ry = [cos(theta/2) -sin(theta/2); sin(theta/2) cos(theta/2)]
        
        U_inj = Ry
        for i in 2:N; U_inj = kron([1 0; 0 1], U_inj); end
        
        # Evolución de rho: rho = U_inj * rho * U_inj†
        rho = U_inj * rho * adjoint(U_inj)
        # 2. EVOLUCIÓN UNITARIA EXACTA: rho = U * rho * U†
        #rho = U_evol * rho * U_evol_adj
        # En lugar de una sola exponencial, iteramos n_substeps
        
        
        #Rk4
        n_substeps = 100     # <--- DEFINIDO: Pasos internos de integración
        dt = T_evol / n_substeps # <--- DEFINIDO: Tamaño del paso temporal
        for _ in 1:n_substeps
            rho = step_rk4_matrix(rho, H_dense, dt)
        end
        #Rk4
        # 3. MEDICIÓN: <Z> = Tr(rho * Z)
        for site in 1:N
            history[k, site] = real(tr(rho * Z_ops[site]))
        end
        
        if k % 10 == 0; print("\rStep $k/100"); end
    end

    # E. GRÁFICA CORREGIDA
    println("\n📊 Generando gráfica...")
    plot_expectation_evolution_easy(1:steps, history, N, inputs)
    println("✅ ¡Listo! Imagen guardada como 'foundation_plot_rho.png'")
end

run_nathan_task()