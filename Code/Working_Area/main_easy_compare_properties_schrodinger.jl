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
gamma = 0.005       
n_substeps = 100
dt = T_evol / n_substeps
Experiment_name = "Schrodinger_Full_Tomography_XYZ_projective_true"
Experiment_path = joinpath(SCRIPT_DIR, "../$Experiment_name")
projective_mode = true
gamma_value = 0.5
if !isdir(Experiment_path); mkpath(Experiment_path); end


# ==============================================================================
# 4. SIMULACIÓN PRINCIPAL
# ==============================================================================
# ==============================================================================
# 4. SIMULACIÓN PRINCIPAL (CORREGIDA)
# ==============================================================================
function run_schrodinger_esp_task()
    println("🚀 Iniciando experimento: $Experiment_name")
    
    # A. Setup
    H_op = build_nathan_all_to_all_XX(N, h_val) # Usa tu función real
    H_dense = operator_to_dense_matrix(H_op, N)
    rho_A = operator_to_dense_matrix(initial_state_all_zeros(N), N)
    rho_B = operator_to_dense_matrix(initial_state_all_ones(N), N)
    Random.seed!(1234); inputs = rand(0:1, steps)

    # B. Observables X, Y, Z
    println("📦 Calculando observables X, Y, Z...")
    obs_matrices = Dict{String, Matrix{ComplexF64}}()
    all_labels = String[]
    #Pauli matrices
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
        # 2-body (Para ESP Full)
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

    # C. Dinámica
    println("🔄 Ejecutando dinámica...")
    for k in 1:steps
        s_k = Float64(inputs[k])
        rz = 1.0 - 2.0 * s_k; rx = 0.0; ry = 0.0 # Inyección Z Ortogonal
        
        # INYECCIÓN (Tu función real)
        rho_A = inject_state_EraseWrite_matrix(rho_A, 0, rz, rx, ry)
        rho_B = inject_state_EraseWrite_matrix(rho_B, 0, rz, rx, ry)

        # EVOLUCIÓN RK4
        for _ in 1:n_substeps
            rho_A = step_rk4_matrix(rho_A, H_dense, dt)
            rho_B = step_rk4_matrix(rho_B, H_dense, dt)
        end
        
        # DEPHASING Z
        rho_A = apply_global_dephasing_schrodinger(rho_A, gamma, "Z")
        rho_B = apply_global_dephasing_schrodinger(rho_B, gamma, "Z")

        # MEDIR
        d_sq = 0.0

        for (lbl, op) in obs_matrices
            # --- NUEVA LLAMADA ---
            # 1. Medimos A
            val_A, rho_A_new = measure_observable_matrix(rho_A, op, projective=projective_mode,gamma_value)
            
            # 2. Medimos B
            val_B, rho_B_new = measure_observable_matrix(rho_B, op, projective=projective_mode,gamma_value)
            ##### Medida sin romper
            #val_A = real(tr(rho_A * op)); val_B = real(tr(rho_B * op))
            # 3. Guardamos datos
            dict_A[lbl][k] = val_A
            dict_B[lbl][k] = val_B
            d_sq += (val_A - val_B)^2
            
            # 4. ACTUALIZAR RHO (Solo si hay backaction)
            # Si projective_mode es true, esto actualiza el estado al colapsado.
            # Si es false, rho_A_new es igual al rho_A viejo, así que no pasa nada.
            if projective_mode
                rho_A = rho_A_new
                rho_B = rho_B_new
            end
        end
        separation_dist[k] = sqrt(d_sq)
        if k%10==0; print("\rStep $k"); end
    end

    # D. PLOTTING (SIN TRUCOS DE RENOMBRADO)
    println("\n📊 Generando gráficas separadas...")
    
    for basis in ["X", "Y", "Z"]
        println("   🎨 Procesando Base $basis...")
        basis_dir = joinpath(Experiment_path, "Basis_$basis")
        if !isdir(basis_dir); mkpath(basis_dir); end
        
        # Filtramos los datos REALES de la base actual
        # (Sin cambiarles el nombre a Z)
        sub_dict_A = filter(p -> contains(p.first, basis), dict_A)
        sub_dict_B = filter(p -> contains(p.first, basis), dict_B)
        
        try
            # Como las funciones ahora son INTELIGENTES (Auto-Detect),
            # entienden que "1X1111" es válido.
            plot_quick_validation_per_qubit(separation_dist, sub_dict_B, inputs, N, basis_dir)
            plot_all_qubits_scatter(sub_dict_B, inputs, N, basis_dir)
            plot_separability_boxplots(sub_dict_A, inputs, N, basis_dir)
            
            # Y tu función Full también (asegúrate de usar la versión corregida que te pasé antes)
            plot_and_save_validation_full(sub_dict_A, sub_dict_B, separation_dist, N, steps, basis_dir)
            
        catch e
            println("Error en plots $basis: $e")
            Base.showerror(stdout, e, catch_backtrace())
        end
    end
    println("\n✅ Experimento finalizado.")
    #SAVE
    save_qrc_results_jld2(N, steps, T_evol, h_val, inputs, dict_A, Experiment_name)
end

# Ejecutar (asegúrate de que las funciones step_rk4 y demás están definidas antes)
run_schrodinger_esp_task()