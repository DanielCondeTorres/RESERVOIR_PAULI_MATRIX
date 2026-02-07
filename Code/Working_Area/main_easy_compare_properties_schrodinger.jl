using LinearAlgebra
using Plots
using Random
using JLD2
using FileIO
using Statistics

# ==============================================================================
# 1. CONFIGURACIÓN DE RUTAS Y CARGA
# ==============================================================================
const SCRIPT_DIR = @__DIR__
function include_rel(path...)
    include(joinpath(SCRIPT_DIR, path...))
end

# Intentamos cargar tus módulos (con fallback si fallan)
try
    include_rel("src/operator_terms/pauli_algebra.jl")
    include_rel("src/operator_terms/hamiltonian.jl")
    include_rel("src/utils/initial_state.jl")
    include_rel("src/utils/measurements.jl")
    include_rel("src/utils/dynamics.jl")
    include_rel("src/auxiliary_scripts/save_files_like_nathan.jl")
    include_rel("src/utils/injection_EraseWrite.jl")
    include_rel("src/visualization/validation_plots.jl")
    include_rel("src/utils/quantum_channels.jl")
    include_rel("src/metrics/separability_metrics.jl")
    include_rel("src/visualization/separability_plots.jl")
    include_rel("src/visualization/separability_and_Echo_property.jl")


catch e
    println("⚠️  Nota: Algunos módulos externos no cargaron. Usando funciones locales.")
end

# ==============================================================================
# 2. PARÁMETROS OPTIMIZADOS
# ==============================================================================
N = 6
steps = 100
T_evol = 1.0        # Tiempo de mezcla
h_val = 1.0
gamma = 0.001        # Dephasing ajustado para ESP
n_substeps = 100
dt = T_evol / n_substeps
Experiment_name = "Schrodinger_EraseWrite_Fixed_Full"
Experiment_path = joinpath(SCRIPT_DIR, "../$Experiment_name")

if !isdir(Experiment_path); mkpath(Experiment_path); end


# ==============================================================================
# 4. SIMULACIÓN PRINCIPAL
# ==============================================================================
function run_schrodinger_esp_task()
    println("🚀 Iniciando experimento: $Experiment_name")
    
    # A. Hamiltoniano
    H_op = build_nathan_all_to_all_XX(N, h_val)
    H_dense = operator_to_dense_matrix(H_op, N)
    U_evol = exp(-im * H_dense * T_evol)
    U_adj = U_evol'

    # B. Estados Iniciales

    rho_A = initial_state_all_zeros(N)
    rho_B = initial_state_all_ones(N)
    rho_A = operator_to_dense_matrix(rho_A, N)
    rho_B = operator_to_dense_matrix(rho_B, N)

    Random.seed!(1234)
    inputs = rand(0:1, steps)

    # C. OBSERVABLES (Z + ZZ) - RESTAURADO!
    println("📦 Calculando operadores observables (Z y ZZ)...")
    obs_matrices = Dict{String, Matrix{ComplexF64}}()
    all_labels = String[]
    sz = [1 0; 0 -1]; id = [1 0; 0 1]
    
    # 1. Z Individuales
    for i in 1:N
        lbl = ["1" for _ in 1:N]; lbl[i] = "Z"; push!(all_labels, join(lbl))
        op = (i==1) ? sz : id; for k in 2:N; op = kron(op, (k==i) ? sz : id); end
        obs_matrices[join(lbl)] = op
    end
    
    # 2. Pares ZZ (NECESARIO PARA TUS PLOTS)
    for i in 1:N, j in (i+1):N
        lbl = ["1" for _ in 1:N]; lbl[i] = "Z"; lbl[j] = "Z"; push!(all_labels, join(lbl))
        op = (i==1) ? sz : id
        for k in 2:N; op = kron(op, (k==i || k==j) ? sz : id); end
        obs_matrices[join(lbl)] = op
    end

    dict_A = Dict(lbl => zeros(Float64, steps) for lbl in all_labels)
    dict_B = Dict(lbl => zeros(Float64, steps) for lbl in all_labels)
    separation_dist = zeros(Float64, steps)

    # D. Bucle Principal
    println("🔄 Ejecutando dinámica...")
    for k in 1:steps
        s_k = Float64(inputs[k])

        # Input 0: |0> (Norte, Z=1) | Input 1: |+> (Ecuador, Z=0)
        theta = (s_k == 0) ? 0.0 : pi/2.0
        rz = cos(theta); rx = sin(theta)
        # 1. Inyección 
        rho_A = inject_state_EraseWrite_matrix(rho_A, 0, rz, rx,0.0)
        rho_B = inject_state_EraseWrite_matrix(rho_B, 0, rz, rx,0.0)

        # 2. Evolución
        #rho_A = U_evol * rho_A * U_adj
        #rho_B = U_evol * rho_B * U_adj
        # 2. EVOLUCIÓN CON RK4
        for _ in 1:n_substeps
            rho_A = step_rk4_matrix(rho_A, H_dense, dt)
            rho_B = step_rk4_matrix(rho_B, H_dense, dt)
        end
        # 3. Dephasing
        rho_A = apply_global_dephasing_schrodinger(rho_A, gamma,"Z")
        rho_B = apply_global_dephasing_schrodinger(rho_B, gamma,"Z")

        # 4. Medir
        dist_sq = 0.0
        for (lbl, op) in obs_matrices
            val_A = real(tr(rho_A * op))
            val_B = real(tr(rho_B * op))
            dict_A[lbl][k] = val_A
            dict_B[lbl][k] = val_B
            dist_sq += (val_A - val_B)^2
        end
        separation_dist[k] = sqrt(dist_sq)
        
        if k % 10 == 0; print("\rStep $k/$steps | Dist: $(round(separation_dist[k], digits=4))"); end
    end

    # E. Gráficas
    println("\n📊 Generando gráficas...")
    plot_quick_validation_per_qubit(separation_dist, dict_B, inputs, N, Experiment_path)
    plot_all_qubits_scatter(dict_B, inputs, N, Experiment_path)
    # 2. Plots Completos (Guardados en disco)
    println("💾 Guardando plots completos en: $Experiment_path")
    try
        plot_and_save_validation_full(dict_A, dict_B, separation_dist, N, steps, Experiment_path)
        plot_separability_boxplots(dict_A, inputs, N, Experiment_path)
    catch e
        println("⚠️ Error en funciones externas de plot: $e")
    end
    
    println("✅ ¡Experimento Completado!")
end

run_schrodinger_esp_task()