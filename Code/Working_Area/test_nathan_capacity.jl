using LinearAlgebra
using Plots
using Random
using JLD2
using FileIO
using Statistics

# ==============================================================================
# 1. CONFIGURACIÓN DE RUTAS Y CARGA DE MÓDULOS
# ==============================================================================
const SCRIPT_DIR = @__DIR__
function include_rel(path...)
    include(joinpath(SCRIPT_DIR, path...))
end

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
    include_rel("src/capacity_training/qrc_training.jl")
    include_rel("src/capacity_training/calculate_stm_capacity.jl")
catch e
    println("⚠️  Nota: Algunos módulos no cargaron. Asegúrate de que las rutas 'src/...' sean correctas.")
end

# ==============================================================================
# 2. FUNCIONES DE DEPHASING (Versión Corregida)
# ==============================================================================

function rotate_density_matrix(rho, axis, n_qubits; inverse=false)
    axis_norm = uppercase(axis)
    if axis_norm == "Z"; return rho; end
    
    gate = if axis_norm == "X"
        1/sqrt(2) * [1.0 1.0; 1.0 -1.0] # Hadamard
    elseif axis_norm == "Y"
        1/sqrt(2) * [1.0 -1.0im; -1.0im 1.0] # Rx(pi/2)
    else
        error("Eje no válido: $axis")
    end
    
    if inverse; gate = gate'; end
    
    U_global = gate
    for i in 2:n_qubits
        U_global = kron(U_global, gate)
    end
    return U_global * rho * U_global'
end

function apply_global_dephasing_matrix(rho::Matrix{ComplexF64}, g::Float64, axis::String="Z")
    if g <= 1e-9; return rho; end
    dim = size(rho, 1)
    n_qubits = Int(log2(dim))
    decay_factor = exp(-(g^2) / 2.0) 
    
    # 1. Rotar al eje de ruido
    rho_rotated = rotate_density_matrix(rho, axis, n_qubits; inverse=false)
    
    # 2. Aplicar decaimiento off-diagonal
    new_rho = copy(rho_rotated)
    for c in 1:dim, r in 1:dim
        if r != c
            n_diff = count_ones((r-1) ⊻ (c-1)) 
            new_rho[r, c] *= (decay_factor ^ n_diff)
        end
    end
    
    # 3. Volver a la base original
    return rotate_density_matrix(new_rho, axis, n_qubits; inverse=true)
end

# ==============================================================================
# 3. PARÁMETROS DEL EXPERIMENTO
# ==============================================================================
N = 6
T_evol = 1.0        
h_val = 1.0
gamma = 0.04      
n_substeps = 100
dt = T_evol / n_substeps
projective_mode = false
gamma_value = 0.8

# Configuración de Archivos
INPUT_FILE = "6_1_2_all_zeros_12345.jld2"
Experiment_name = "Schrodinger_Replay_s_vec"
Experiment_path = joinpath(SCRIPT_DIR, "../$Experiment_name")
if !isdir(Experiment_path); mkpath(Experiment_path); end

# Parámetros STM
WASHOUT = 50      
MAX_DELAY = 15     

# ==============================================================================
# 4. SIMULACIÓN PRINCIPAL
# ==============================================================================
function run_schrodinger_esp_task()
    println("🚀 Iniciando experimento: $Experiment_name")
    
    # --- CARGA DE INPUTS (Buscando s_vec) ---
    input_file_path = isfile(INPUT_FILE) ? INPUT_FILE : joinpath(SCRIPT_DIR, "../Input_Data", INPUT_FILE)
    if !isfile(input_file_path); error("❌ No se encuentra el archivo de inputs: $input_file_path"); end
    
    println("📂 Cargando inputs desde: $input_file_path")
    data_loaded = load(input_file_path)
    
    # Extraemos s_vec del diccionario meta_data_dict
    inputs = if haskey(data_loaded, "meta_data_dict") && haskey(data_loaded["meta_data_dict"], "s_vec")
        data_loaded["meta_data_dict"]["s_vec"]
    else
        error("❌ No se encontró 's_vec' dentro de 'meta_data_dict'. Claves: $(keys(data_loaded))")
    end
    
    steps = length(inputs)
    println("✅ Configuración: $N qubits, $steps pasos detectados.")

    # --- SETUP SISTEMA ---
    H_op = build_nathan_all_to_all_XX(N, h_val) 
    H_dense = operator_to_dense_matrix(H_op, N)
    rho_A = operator_to_dense_matrix(initial_state_all_zeros(N), N)
    rho_B = operator_to_dense_matrix(initial_state_all_ones(N), N)

    # --- DEFINICIÓN DE OBSERVABLES ---
    obs_matrices = Dict{String, Matrix{ComplexF64}}()
    all_labels = String[]
    sx=[0. 1.; 1. 0.]; sy=[0. -im; im 0.]; sz=[1. 0.; 0. -1.]; id=[1. 0.; 0. 1.]
    
    for (b_char, b_op) in zip(['X','Y','Z'], [sx, sy, sz])
        b_str = string(b_char)
        for i in 1:N
            lbl = ["1" for _ in 1:N]; lbl[i] = b_str; s_lbl = join(lbl)
            push!(all_labels, s_lbl)
            op = (i==1) ? b_op : id; for k in 2:N; op = kron(op, (k==i) ? b_op : id); end
            obs_matrices[s_lbl] = op
        end
        for i in 1:N, j in (i+1):N
            lbl = ["1" for _ in 1:N]; lbl[i]=b_str; lbl[j]=b_str; s_lbl = join(lbl)
            push!(all_labels, s_lbl)
            op = (i==1) ? b_op : id; for k in 2:N; op = kron(op, (k==i||k==j) ? b_op : id); end
            obs_matrices[s_lbl] = op
        end
    end

    dict_A = Dict(lbl => zeros(Float64, steps) for lbl in all_labels)
    dict_B = Dict(lbl => zeros(Float64, steps) for lbl in all_labels)
    separation_dist = zeros(Float64, steps)

    # --- BUCLE DE DINÁMICA ---
    println("🔄 Ejecutando dinámica...")
    for k in 1:steps
        s_k = Float64(inputs[k])
        rz = 1.0 - 2.0 * s_k; rx = 0.0; ry = 0.0 
        
        rho_A = inject_state_EraseWrite_matrix(rho_A, 0, rz, rx, ry)
        rho_B = inject_state_EraseWrite_matrix(rho_B, 0, rz, rx, ry)

        for _ in 1:n_substeps
            rho_A = step_rk4_matrix(rho_A, H_dense, dt)
            rho_B = step_rk4_matrix(rho_B, H_dense, dt)
        end
        
        rho_A = apply_global_dephasing_matrix(rho_A, gamma, "Z")
        rho_B = apply_global_dephasing_matrix(rho_B, gamma, "Z")

        d_sq = 0.0
        for (lbl, op) in obs_matrices
            val_A, rho_A_new = measure_observable_matrix(rho_A, op, projective=projective_mode, gamma_value)
            val_B, rho_B_new = measure_observable_matrix(rho_B, op, projective=projective_mode, gamma_value)
            
            dict_A[lbl][k] = val_A
            dict_B[lbl][k] = val_B
            d_sq += (val_A - val_B)^2
            
            if projective_mode
                rho_A, rho_B = rho_A_new, rho_B_new
            end
        end
        separation_dist[k] = sqrt(d_sq)
        if k%10==0; print("\rStep $k / $steps"); end
    end

    # --- POST-PROCESADO Y STM ---
    println("\n🧠 Calculando Capacidad de Memoria...")
    # La capacidad usa x_k * W para predecir s_{k-tau}
    capacities, total_stm = calculate_stm_capacity(dict_B, inputs, MAX_DELAY, WASHOUT)
    plot_memory_capacity(capacities, MAX_DELAY, total_stm, Experiment_path)

    println("📊 Generando Boxplots y Validaciones...")
    for basis in ["X", "Y", "Z"]
        basis_dir = joinpath(Experiment_path, "Basis_$basis"); mkpath(basis_dir)
        sub_A = filter(p -> contains(p.first, basis), dict_A)
        sub_B = filter(p -> contains(p.first, basis), dict_B)
        plot_separability_boxplots(filter(p -> contains(p.first, basis), dict_A), inputs, N, basis_dir)
        plot_quick_validation_per_qubit(separation_dist, sub_B, inputs, N, basis_dir)
        plot_and_save_validation_full(sub_A, sub_B, separation_dist, N, steps, basis_dir)

    end

    save_qrc_results_jld2(N, steps, T_evol, h_val, inputs, dict_B, Experiment_name)
    println("\n✅ Finalizado. Total STM: $total_stm")
end

run_schrodinger_esp_task()