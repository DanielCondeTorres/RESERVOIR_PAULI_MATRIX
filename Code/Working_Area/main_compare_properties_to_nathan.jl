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
    println("⚠️  Nota: Algunos módulos externos no cargaron. Usando funciones locales.")
end

# ==============================================================================
# 2. PARÁMETROS
# ==============================================================================
N = 6
steps = 100
T_evol = 1.0        
h_val = 1.0
gamma = 0.04      
n_substeps = 100
dt = T_evol / n_substeps
Experiment_name = "Schrodinger_Full_Tomography_XYZ_projective_false_continuo_2"
Experiment_path = joinpath(SCRIPT_DIR, "../$Experiment_name")
projective_mode = false
gamma_value = 0.8
if !isdir(Experiment_path); mkpath(Experiment_path); end

# Parámetros para la Capacidad
WASHOUT = 2      
MAX_DELAY = 10     
TRAIN_RATIO = 0.9  

# ==============================================================================
# 4. SIMULACIÓN PRINCIPAL
# ==============================================================================
function run_schrodinger_esp_task()
    println("🚀 Iniciando experimento: $Experiment_name")
    
    # A. Setup
    H_op = build_nathan_all_to_all_XX(N, h_val)
    H_dense = operator_to_dense_matrix(H_op, N)
    U_evol = exp(-im * H_dense * T_evol)
    U_evol_adj = adjoint(U_evol)
    
    rho_A = operator_to_dense_matrix(initial_state_all_zeros(N), N)
    rho_B = operator_to_dense_matrix(initial_state_all_ones(N), N)
    rho_C = operator_to_dense_matrix(initial_state_all_zeros(N), N)
    
    Random.seed!(1234)
    inputs = rand(steps)
    
    inputs_C = copy(inputs)
    step_perturbation = 30
    if steps > step_perturbation
        inputs_C[step_perturbation+1:end] = rand(steps - step_perturbation)
    end

    # B. Observables X, Y, Z (1, 2 y 3 cuerpos)
    println("📦 Calculando observables X, Y, Z (incluyendo tríadas)...")
    obs_matrices = Dict{String, Matrix{ComplexF64}}()
    all_labels = String[]
    sx=[0. 1.; 1. 0.]; sy=[0. -im; im 0.]; sz=[1. 0.; 0. -1.]; id=[1. 0.; 0. 1.]
    
    for (b_char, b_op) in zip(['X','Y','Z'], [sx, sy, sz])
        b_str = string(b_char)
        
        # 1-body
        for i in 1:N
            lbl = ["1" for _ in 1:N]; lbl[i] = b_str; s_lbl = join(lbl)
            push!(all_labels, s_lbl)
            op = (i==1) ? b_op : id; for k in 2:N; op = kron(op, (k==i) ? b_op : id); end
            obs_matrices[s_lbl] = op
        end
        
        # 2-body
        for i in 1:N, j in (i+1):N
            lbl = ["1" for _ in 1:N]; lbl[i]=b_str; lbl[j]=b_str; s_lbl = join(lbl)
            push!(all_labels, s_lbl)
            op = (i==1) ? b_op : id; for k in 2:N; op = kron(op, (k==i||k==j) ? b_op : id); end
            obs_matrices[s_lbl] = op
        end

        # 3-body (TRÍADAS) <--- NUEVO BLOQUE
        for i in 1:N, j in (i+1):N, k in (j+1):N
            lbl = ["1" for _ in 1:N]
            lbl[i] = b_str; lbl[j] = b_str; lbl[k] = b_str
            s_lbl = join(lbl)
            push!(all_labels, s_lbl)
            
            # Construcción del operador: Kronecker product
            op = (1==i || 1==j || 1==k) ? b_op : id
            for m in 2:N
                op = kron(op, (m==i || m==j || m==k) ? b_op : id)
            end
            obs_matrices[s_lbl] = op
        end
    end

    # Inicialización de diccionarios
    dict_A = Dict(lbl => zeros(Float64, steps) for lbl in all_labels)
    dict_B = Dict(lbl => zeros(Float64, steps) for lbl in all_labels)
    dict_C = Dict(lbl => zeros(Float64, steps) for lbl in all_labels)
    
    separation_dist = zeros(Float64, steps)      
    separation_dist_AC = zeros(Float64, steps)   

    # C. Dinámica
    println("🔄 Ejecutando dinámica...")
    for k in 1:steps
        s_k = Float64(inputs[k])
        rz = 1.0 - 2.0 * s_k
        
        s_k_C = Float64(inputs_C[k])
        rz_C = 1.0 - 2.0 * s_k_C
        
        rx = 0.0; ry = 0.0 
        
        rho_A = inject_state_EraseWrite_matrix(rho_A, 0, rz, rx, ry)
        rho_B = inject_state_EraseWrite_matrix(rho_B, 0, rz, rx, ry)
        rho_C = inject_state_EraseWrite_matrix(rho_C, 0, rz_C, rx, ry)

        for _ in 1:n_substeps
            rho_A = U_evol * rho_A * U_evol_adj
            rho_B = U_evol * rho_B * U_evol_adj
            rho_C = U_evol * rho_C * U_evol_adj
        end
        
        d_sq = 0.0
        d_sq_AC = 0.0

        for (lbl, op) in obs_matrices
            val_A, _ = measure_observable_matrix(rho_A, op, projective=projective_mode, gamma_value)
            val_B, _ = measure_observable_matrix(rho_B, op, projective=projective_mode, gamma_value)
            val_C, _ = measure_observable_matrix(rho_C, op, projective=projective_mode, gamma_value)
            
            dict_A[lbl][k] = val_A
            dict_B[lbl][k] = val_B
            dict_C[lbl][k] = val_C

            # Sumamos a la distancia (usando todos los observables, incluidos tríadas)
            d_sq += (val_A - val_B)^2       
            d_sq_AC += (val_A - val_C)^2    
        end
        separation_dist[k] = sqrt(d_sq)
        separation_dist_AC[k] = sqrt(d_sq_AC)
        
        if k%10==0; print("\rStep $k"); end
    end

    # D. PLOTTING
    println("\n📊 Generando gráficas separadas...")
    
    for basis in ["X", "Y", "Z"]
        println("   🎨 Procesando Base $basis...")
        basis_dir = joinpath(Experiment_path, "Basis_$basis")
        if !isdir(basis_dir); mkpath(basis_dir); end
        
        sub_dict_A = filter(p -> contains(p.first, basis), dict_A)
        sub_dict_B = filter(p -> contains(p.first, basis), dict_B)
        sub_dict_C = filter(p -> contains(p.first, basis), dict_C)

        try
            # 1. Plots originales (A vs B)
            # Nota: plot_quick_validation y scatter podrían saturarse con las tríadas, 
            # pero funcionarán si están diseñadas genéricamente.
            
            # AQUI: Llamada normal A vs B con nueva función de plot
            plot_and_save_validation_full(
                sub_dict_A, sub_dict_B, separation_dist, N, steps, basis_dir;
                label_A="Tray. A", 
                label_B="Tray. B"
            )
            
            # 2. NUEVO PLOT (A vs C)
            dir_AC = joinpath(basis_dir, "Compare_A_C")
            mkpath(dir_AC)
            println("      Guardando comparación A vs C en: $dir_AC")
            
            plot_and_save_validation_full(
                sub_dict_A, sub_dict_C, separation_dist_AC, N, steps, dir_AC;
                label_A="Tray. A (Orig)", 
                label_B="Tray. C (Perturb)"
            )
            
        catch e
            println("Error en plots $basis: $e")
            Base.showerror(stdout, e, catch_backtrace())
        end
    end
    println("\n✅ Experimento finalizado.")
    
    save_qrc_results_jld2(N, steps, T_evol, h_val, inputs, dict_A, Experiment_name)

    capacities, total_stm = calculate_stm_capacity(dict_B, inputs, MAX_DELAY, 50)
    plot_memory_capacity(capacities, MAX_DELAY, total_stm, Experiment_path)
end

run_schrodinger_esp_task()