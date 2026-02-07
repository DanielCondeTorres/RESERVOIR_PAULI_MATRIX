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
gamma = 0.05        # Dephasing ajustado para ESP
Experiment_name = "Schrodinger_EraseWrite_Fixed_Full"
Experiment_path = joinpath(SCRIPT_DIR, "../$Experiment_name")

if !isdir(Experiment_path); mkpath(Experiment_path); end

# ==============================================================================
# 3. FUNCIONES AUXILIARES LOCALES (Salvavidas)
# ==============================================================================

# Generador Hamiltoniano Denso
function build_hamiltonian_dense_local(N, h)
    sx = [0 1; 1 0]; sz = [1 0; 0 -1]; id = [1 0; 0 1]
    H = zeros(ComplexF64, 2^N, 2^N)
    # Campo Z
    for i in 1:N
        op = (i==1) ? sz : id; for j in 2:N; op = kron(op, (j==i) ? sz : id); end
        H .+= h * op
    end
    # XX interactions
    Random.seed!(42)
    J = randn(N,N) ./ sqrt(N-1)
    for i in 1:N, j in (i+1):N
        op = (i==1) ? sx : id; for k in 2:N; op = kron(op, (k==i || k==j) ? sx : id); end
        H .+= J[i,j] * op
    end
    return H
end

# Dephasing Z Denso (Optimizado)
function apply_dephasing_local(rho::Matrix{ComplexF64}, g::Float64)
    if g <= 1e-9; return rho; end
    dim = size(rho, 1)
    decay = exp(-g)
    new_rho = copy(rho)
    for c in 1:dim
        for r in 1:dim
            if r == c; continue; end
            # Hamming distance trick
            n_diff = count_ones((r-1) ⊻ (c-1))
            if n_diff > 0; new_rho[r,c] *= (decay^n_diff); end
        end
    end
    return new_rho
end

# Inyección Erase & Write (Matriz)
function inject_ew_local(rho::Matrix{ComplexF64}, qubit_idx::Int, rz::Float64, rx::Float64)
    dim = size(rho, 1); Nq = Int(log2(dim))
    I2 = [1.0+0im 0.0; 0.0 1.0]; P0 = [1.0+0im 0.0; 0.0 0.0]; Flip = [0.0+0im 1.0; 0.0 0.0]
    
    k = qubit_idx + 1
    M0 = (k==1) ? P0 : I2; M1 = (k==1) ? Flip : I2
    for i in 2:Nq
        M0 = kron(M0, (i==k) ? P0 : I2)
        M1 = kron(M1, (i==k) ? Flip : I2)
    end
    rho_reset = M0 * rho * M0' + M1 * rho * M1'
    
    # Rotación
    norm = sqrt(rx^2 + rz^2)
    if norm < 1e-9; return rho_reset; end
    theta = acos(rz/norm) # Asumiendo ry=0
    Ry = [cos(theta/2) -sin(theta/2); sin(theta/2) cos(theta/2)]
    
    U_op = (k==1) ? Ry : I2
    for i in 2:Nq; U_op = kron(U_op, (i==k) ? Ry : I2); end
    
    return U_op * rho_reset * U_op'
end

# ==============================================================================
# 4. SIMULACIÓN PRINCIPAL
# ==============================================================================
function run_schrodinger_esp_task()
    println("🚀 Iniciando experimento: $Experiment_name")
    
    # A. Hamiltoniano
    H_dense = build_hamiltonian_dense_local(N, h_val)
    U_evol = exp(-im * H_dense * T_evol)
    U_adj = U_evol'

    # B. Estados Iniciales
    v0 = [1.0+0im; 0.0]; psi0 = v0; for i in 2:N; psi0 = kron(psi0, v0); end
    rho_A = psi0 * psi0'
    
    v1 = [0.0+0im; 1.0]; psi1 = v1; for i in 2:N; psi1 = kron(psi1, v1); end
    rho_B = psi1 * psi1'

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
        rho_A = inject_ew_local(rho_A, 0, rz, rx)
        rho_B = inject_ew_local(rho_B, 0, rz, rx)

        # 2. Evolución
        rho_A = U_evol * rho_A * U_adj
        rho_B = U_evol * rho_B * U_adj

        # 3. Dephasing
        rho_A = apply_dephasing_local(rho_A, gamma)
        rho_B = apply_dephasing_local(rho_B, gamma)

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
    plot_quick_validation_per_qubit(separation_dist, dict_A, inputs, N, Experiment_path)
    plot_all_qubits_scatter(dict_A, inputs, N, Experiment_path)
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