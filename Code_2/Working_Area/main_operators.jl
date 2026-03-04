#LOAD MODULES
using LinearAlgebra
using Plots
using Random
using JLD2
using FileIO
using Statistics
using JSON
using Plots
const SCRIPT_DIR = @__DIR__
function include_rel(path...)
    include(joinpath(SCRIPT_DIR, path...))
end
# Cargamos tus archivos
include_rel("src/utils/operators_operations.jl")
include_rel("src/utils/truncate_operators.jl")

include_rel("src/dynamics/runge_kutta_4.jl")   
include_rel("src/dynamics/initial_state.jl")  
include_rel("src/dynamics/hamiltonian.jl")    
include_rel("src/dynamics/global_dephasing_measurement.jl")    
include_rel("src/dynamics/eraser_and_write.jl")


include_rel("src/saver/save_results.jl")    
include_rel("src/saver/save_plots.jl")        



include_rel("src/visualization/visualization_comparision.jl")
include_rel("src/utils/capacity.jl")












# ==============================================================================
# 1. PARAMETERS
# ==============================================================================
# Load from files
#INPUT_FILE = "6_3_2_all_zeros_12345.jld2"
# Parameters to calculate the capacity
# --- CARGA DE INPUTS (Buscando s_vec) ---
#input_file_path = isfile(INPUT_FILE) ? INPUT_FILE : joinpath(SCRIPT_DIR, "../Input_Data/", INPUT_FILE) #"../Input_Data/" is the path to INPUT_FILE
#if !isfile(input_file_path); error("❌ No se encuentra el archivo de inputs: $input_file_path"); end
#println(" Cargando inputs desde: $input_file_path")
#data_loaded = load(input_file_path)
#meta = data_loaded["meta_data_dict"]
#nathan_results = data_loaded["expect_dict"]
# End load from files (can be commented if you want to use random inputs)


#NEW INPUT FILES USING JSON (for better readability), BECAUSE PREVIOUS PROBLEMS

json_raw = "/Users/danielcondetorres/Desktop/IBM_EXAMENES/TFM_NATHAN/TFM_VERSION_2/Code/Input_Data/results.json"
datos = JSON.parsefile(json_raw)

# 2. Navegar hasta la meta-data
# Obtenemos la primera llave (ej: "g=0.3_seed=1...")
llave_raiz = collect(keys(datos["data"]))[1]
meta = datos["data"][llave_raiz]["meta"]

#Save files:
#Save files:
Experiment_name = "EVALUACION_RK4_SIN_ERA_WRITE_SIN_DEPHASING_SIN_TRUNCATE" # Where i save the files
# Usamos abspath para forzar la ruta completa desde la raíz de tu Mac
Experiment_path = abspath(joinpath(SCRIPT_DIR, "../$Experiment_name"))
projective_mode = false # Not using it
if !isdir(Experiment_path); mkpath(Experiment_path); end

println("📁 TODAS LAS GRÁFICAS Y DATOS IRÁN A: $Experiment_path")
# --- SETUP SISTEMA (Usando parámetros del archivo) ---
T_evol = 10.0        
Js= 1.
N = 6#Int(meta["N"])
h_val = Float64(meta["h"])
gamma = Float64(meta["g"]) 
J_vec_file = nothing  # Vector{Float64}(meta["Jvec"])
Random.seed!(1234)
steps = 100
inputs = Float64.(rand(0:1, steps))
#inputs=rand(steps)
n_substeps = 1000  
#inputs = Vector{Float64}(meta["s_vec"]) 








# 3. ASIGNACIÓN AUTOMÁTICA DE TUS VARIABLES
N = 6#Int(meta["N"])               # Ahora será 4 automáticamente
h_val = Float64(meta["h"])       # 10.0
gamma = Float64(meta["g"])       # 0.3
T_evol = Float64(meta["Δt"])     # 10.0
#steps = Int(meta["num_steps"])   # 10
#inputs = Float64.(meta["s_vec"]) # [1.0, 1.0, 1.0, 1.0, 0.0, ...]
steps = 100
inputs=rand(steps)
# Para el J_vec_file (los acoplamientos):
# Como Nathan guardó un Jdict, extraemos solo los valores numéricos
J_dict = meta["Jdict"]
J_vec_file = Float64.(collect(values(J_dict))) # Esto te da [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]



GATE_WHERE_DEPHASING_IS_APPLIED="Z"
println("\n LOADED PARAMETERS:")
println("-"^40)
println("N (Qubits):      $N")
println("Inputs:      $inputs")
println("h (Transversal): $h_val")
println("Gamma (g):       $gamma")
println("Steps (Pasos):   $steps")
#println("dt (Paso temp):  $dt")
# Verificamos si J_vec_file existe antes de pedirle la longitud
longitud_j = isnothing(J_vec_file) ? "0 (None/Random)" : length(J_vec_file)
println("Longitud Jvec:   $longitud_j")
println("Longitud s_vec:  $(length(inputs))")


# ==============================================================================
# 1. FUNCTIONS (EXTERNAL TO THE MAIN)
# ==============================================================================


#PLOTTING ECHO AND SEPARABILITY



function run_task_pauli()
    println("🔮 Iniciando Experimento Pauli-Only: $Experiment_name")
    
    # --- SETUP SISTEMA ---
    # Hamiltoniano nativo en Paulis
    H_op = build_nathan_all_to_all_XX(N, h_val, J_vec_file) 
    
    # Estados iniciales en base Pauli
    rho_A = initial_state_all_zeros(N)
    rho_B = initial_state_all_ones(N)
    rho_C = initial_state_all_zeros(N) # Para chequear separabilidad

    # --- PERTURBACIÓN TRAYECTORIA C ---
    inputs_C = copy(inputs)
    step_perturbation = 50
    if steps > step_perturbation
        # A partir del paso 30, la entrada C diverge de la A
        inputs_C[step_perturbation+1:end] = rand(steps - step_perturbation)
    end

    # --- PREPARACIÓN DE MEDIDAS ---
    all_labels = generate_obs_labels(N)
    dict_A = Dict(lbl => zeros(Float64, steps) for lbl in all_labels)
    dict_B = Dict(lbl => zeros(Float64, steps) for lbl in all_labels)
    dict_C = Dict(lbl => zeros(Float64, steps) for lbl in all_labels)

    # Sub-pasos para estabilidad unitaria en el integrador RK4
    
    dt_sub = T_evol / n_substeps 

    # ==========================================================================
    # 2. BUCLE TEMPORAL (Dinámica del Reservorio)
    # ==========================================================================
    println("🔄 Evolucionando con step_rk4 (Pauli Strings)...")
    for k in 1:steps
        # Encoding: Mapeo del input a la esfera de Bloch (Eje Z)
        u_A = inputs[k];   rz_A = 1.0 - 2.0 * u_A
        u_C = inputs_C[k]; rz_C = 1.0 - 2.0 * u_C
        rx, ry = 0.0, 0.0

        # A. INYECCIÓN: Erase & Write (Reset local en Qubit 1)
        rho_A = inject_state_EraseWrite_pauli(rho_A, 0, rz_A, rx, ry)
        rho_B = inject_state_EraseWrite_pauli(rho_B, 0, rz_A, rx, ry)
        rho_C = inject_state_EraseWrite_pauli(rho_C, 0, rz_C, rx, ry)

        # B. EVOLUCIÓN UNITARIA: Usando tu step_rk4
        for _ in 1:n_substeps
            rho_A = step_rk4(rho_A, H_op, dt_sub)
            rho_B = step_rk4(rho_B, H_op, dt_sub)
            rho_C = step_rk4(rho_C, H_op, dt_sub)

            MAX_TERMS = 1000 # El valor M que decidas (ej. 1000, 5000...)
            rho_A = truncate_operator(rho_A, MAX_TERMS)
            rho_B = truncate_operator(rho_B, MAX_TERMS)
            rho_C = truncate_operator(rho_C, MAX_TERMS)
        end

        # C. Dephasing Global
        rho_A = apply_global_dephasing(rho_A, gamma, GATE_WHERE_DEPHASING_IS_APPLIED)
        rho_B = apply_global_dephasing(rho_B, gamma, GATE_WHERE_DEPHASING_IS_APPLIED)
        rho_C = apply_global_dephasing(rho_C, gamma, GATE_WHERE_DEPHASING_IS_APPLIED)

        # D. MEDIDA: <P> = c_k * 2^N
        for lbl in all_labels
            p_str = label_to_pauli(lbl, N)
            dict_A[lbl][k] = real(get(rho_A, p_str, 0.0im) * (2^N))
            dict_B[lbl][k] = real(get(rho_B, p_str, 0.0im) * (2^N))
            dict_C[lbl][k] = real(get(rho_C, p_str, 0.0im) * (2^N))
        end
        
        if k % 10 == 0; print("\rProgreso: $k/$steps"); end
    end


# ==========================================================================
    # 3. POST-PROCESADO Y PLOTTING
    # ==========================================================================
    println("\n🎨 Generando gráficas de validación...")
    for basis in ["X", "Y", "Z"]
        basis_dir = joinpath(Experiment_path, "Basis_$basis")
        if !isdir(basis_dir); mkpath(basis_dir); end
        
        sub_dict_A = filter(p -> contains(p.first, basis), dict_A)
        sub_dict_B = filter(p -> contains(p.first, basis), dict_B)
        sub_dict_C = filter(p -> contains(p.first, basis), dict_C)

        # Validación ESP (Echo State Property): A vs B
        dir_AB = joinpath(basis_dir, "Check_ECHO_A_B")
        plot_and_save_validation_full(sub_dict_A, sub_dict_B, N, steps, dir_AB; 
                                      label_A="Tray. A", label_B="Tray. B")

        # Validación de Separabilidad: A vs C
        dir_AC = joinpath(basis_dir, "Check_separability_A_C")
        plot_and_save_validation_full(sub_dict_A, sub_dict_C, N, steps, dir_AC; 
                                      label_A="Original (A)", label_B="Perturbada (C)")
    end

    # ¡¡¡¡AQUÍ ESTABA EL ERROR!!!! TIENE QUE SER Experiment_path
    save_qrc_results_jld2(N, steps, T_evol, h_val, inputs, dict_A, Experiment_path)
    println("\n✅ Benchmarking completado correctamente.")



    println("\n📊 Calculando Short-Term Memory Capacity...")
    all_labels_A = collect(keys(dict_A))
    states_matrix_A = zeros(Float64, steps, length(all_labels_A))
    for (j, lbl) in enumerate(all_labels_A)
        states_matrix_A[:, j] .= dict_A[lbl]
    end
    # Using tau_max = 10, washout = 10 (discard first 10 steps)
    capacity_A = calculate_stm_capacity(states_matrix_A, inputs, 10; washout=10)
    # Plot capacity
    plot_stm_capacity(capacity_A, joinpath(Experiment_path, "STM_Capacity.png"))
        
end
plt =  guardar_graficas_individuales(json_raw,Experiment_path)
display(plt)
run_task_pauli()
