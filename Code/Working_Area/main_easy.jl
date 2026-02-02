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

# ==============================================================================
# 2. PARÁMETROS
# ==============================================================================
N = 6
steps = 100
T_evol = 10.0
h_val = 1.0
Experiment_name = "../Nathan_Foundation_AllToAll_FullZ"

# Parámetros para RK4
n_substeps = 100
dt = T_evol / n_substeps

# Cargamos tus archivos
include_rel("src/operator_terms/pauli_algebra.jl")      
include_rel("src/operator_terms/hamiltonian.jl")        
include_rel("src/utils/initial_state.jl") 
include_rel("src/utils/measurements.jl")
include_rel("src/utils/dynamics.jl")
include_rel("src/visualization/expectation_xj_vs_step.jl") 
include_rel("src/auxiliary_scripts/save_files_like_nathan.jl")   



# ==============================================================================
# FUNCIONES AUXILIARES
# ==============================================================================
function dense_to_operator(rho_mat::Matrix{ComplexF64}, n_qubits::Int)
    rho_op = Operator()
    dim = 2^n_qubits
    scale = 1.0 / dim
    for site in 0:(n_qubits-1)
        z_mask = 1 << site
        p = PauliString(0, z_mask)
        P_mat = [1.0;;]
        for i in 0:(n_qubits-1)
            gate = (i == site) ? [1 0; 0 -1] : [1 0; 0 1]
            P_mat = kron(gate, P_mat)
        end 
        rho_op[p] = tr(rho_mat * P_mat) * scale
    end
    return rho_op
end

# ==============================================================================
# 3. SIMULACIÓN EXACTA (RK4 + FULL MEASUREMENTS)
# ==============================================================================
function run_nathan_task()
    println(" Iniciando experimento: $Experiment_name")
    println("Configuración: N=$N, h=$h_val, steps=$steps")

    # A. Construir Hamiltoniano
    H_op = build_nathan_all_to_all_XX(N, h_val)
    H_dense = operator_to_dense_matrix(H_op, N)
   
    
    
    U_evol = exp(-im * H_dense * T_evol)
    U_evol_adj = adjoint(U_evol) # U† para la evolución de rho

    # B. Preparar Estado Inicial e Inputs
    println(" Creando estado inicial |00...0>...")
    rho_op = initial_state_all_zeros(N)
    rho = operator_to_dense_matrix(rho_op, N)
    
    Random.seed!(1234) # Semilla para que los inputs sean reproducibles
    inputs = rand(0:1, steps)

    # C. Pre-calcular TODOS los Observables (Z y ZZ)
    println(" Pre-calculando observables Z_i y correlaciones Z_i Z_j...")
    obs_matrices = Dict{String, Matrix{ComplexF64}}()
    
    # 1. Individuales: "Z11111", "1Z1111"...
    for i in 1:N
        label = ["1" for _ in 1:N]; label[i] = "Z"
        p_string = create_pauli_observable("Z", [i])
        obs_matrices[join(label)] = operator_to_dense_matrix(Operator(p_string => 1.0+0.0im), N)
    end
    
    # 2. Pares ZZ: "ZZ1111", "Z1Z111"...
    for i in 1:N, j in (i+1):N
        label = ["1" for _ in 1:N]; label[i] = "Z"; label[j] = "Z"
        p_string = create_pauli_observable("Z", [i, j])
        obs_matrices[join(label)] = operator_to_dense_matrix(Operator(p_string => 1.0+0.0im), N)
    end

    # Diccionario para guardar toda la historia
    expect_dict = Dict{String, Vector{Float64}}()
    for label in keys(obs_matrices)
        expect_dict[label] = zeros(Float64, steps)
    end

    # D. Bucle de Reservoir
    println(" Ejecutando pasos temporales con RK4...")
    for k in 1:steps
        s_k = Float64(inputs[k])

        # 1. INYECT STATE (Encoding Ry)
        theta = 2.0 * asin(sqrt(s_k))
        Ry = [cos(theta/2) -sin(theta/2); sin(theta/2) cos(theta/2)]
        U_inj = Ry
        for i in 2:N; U_inj = kron([1 0; 0 1], U_inj); end
        
        rho = U_inj * rho * adjoint(U_inj)
        rho = U_evol * rho * U_evol_adj
        # 2. EVOLUCIÓN CON RK4
        #for _ in 1:n_substeps
        #    rho = step_rk4_matrix(rho, H_dense, dt)
        #end

        # 3. MEASURE: Medimos todos los operadores del diccionario
        for (label, matrix_op) in obs_matrices
            expect_dict[label][k] = real(tr(rho * matrix_op))
        end
        
        if k % 10 == 0; print("\rStep $k/$steps"); end
    end

    # ==========================================================================
    # E. GUARDADO Y GRÁFICA
    # ==========================================================================
    # 1. Guardar JLD2 usando tu función modular
    ruta_experimento=save_qrc_results_jld2(N, steps, T_evol, h_val, inputs, expect_dict, Experiment_name)
    
    # 2. Preparar matriz de historia solo para las Z individuales para el plot
    history_single_z = zeros(Float64, steps, N)
    for i in 1:N
        label = ["1" for _ in 1:N]; label[i] = "Z"
        history_single_z[:, i] = expect_dict[join(label)]
    end

    println("\n Generando gráfica...")
    # Llamamos a tu función de plot. Asegúrate de que use Experiment_name para el título/archivo
    plot_expectation_evolution_easy(1:steps, history_single_z, N, ruta_experimento, inputs)
    plot_zz_correlations(1:steps, expect_dict, ruta_experimento, inputs)
    println("¡Listo! Experimento '$Experiment_name' finalizado.")
end
# --- AL FINAL DEL ARCHIVO ---
run_nathan_task()
