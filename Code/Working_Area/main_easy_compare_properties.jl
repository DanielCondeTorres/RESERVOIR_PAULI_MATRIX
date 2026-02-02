using LinearAlgebra
using Plots
using Random
using JLD2
using FileIO
using Dates 

# ==============================================================================
# 1. CONFIGURACIÓN DE RUTAS Y AUXILIARES
# ==============================================================================
const SCRIPT_DIR = @__DIR__
function include_rel(path...)
    include(joinpath(SCRIPT_DIR, path...))
end

# ==============================================================================
# 2. CARGA DE MÓDULOS
# ==============================================================================
include_rel("src/operator_terms/pauli_algebra.jl")      
include_rel("src/operator_terms/hamiltonian.jl")        
include_rel("src/utils/initial_state.jl") 
include_rel("src/utils/measurements.jl")
include_rel("src/utils/dynamics.jl")
include_rel("src/auxiliary_scripts/save_files_like_nathan.jl")   
include_rel("src/utils/injection_EraseWrite.jl")

# NUEVO: Cargamos la función de visualización
include_rel("src/visualization/validation_plots.jl") 
# ==============================================================================
# 3. FUNCIÓN PRINCIPAL DE SIMULACIÓN Y VALIDACIÓN
# ==============================================================================
function run_nathan_validation_task()
    # --- Parámetros ---
    N = 6
    steps = 200
    T_evol = 10.0
    h_val = 10.
    n_substeps = 10000
    dt = T_evol / n_substeps
    g=0.05  # Tasa de dephasing
    # --- RUTA DE GUARDADO (Corregida para que use SCRIPT_DIR) ---
    folder_name = "../RESULTADOS_ESP_SEPARACION"
    Experiment_path = joinpath(SCRIPT_DIR, folder_name)
    
    println("🚀 Iniciando Validación QRC: ESP y Separación")
    println("📂 Creando directorio en: $Experiment_path")
    
    if !isdir(Experiment_path); mkdir(Experiment_path); end

    H_op = build_nathan_all_to_all_XX(N, h_val) 
    
    # Pre-calculamos etiquetas
    z_labels = [join([k == i ? "Z" : "1" for k in 1:N]) for i in 1:N]
    label_to_obj = Dict(lbl => string_to_pauli(lbl) for lbl in z_labels)
    
    # --- SETUP PARA TEST ESP ---
    rho_A = initial_state_all_zeros(N) 
    
    # Rho B diferente: Usamos tu función o rotación manual si no existe
    # rho_B = initial_state_x_first(N) <--- Úsala si la tienes definida
    # Si no, rotación manual para asegurar diferencia:
    rho_B = initial_state_all_zeros(N)
    rho_B = initial_state_all_ones(N)
    
    Random.seed!(1234)
    inputs = rand(0:1, steps)

    hist_A = zeros(steps, N)
    hist_B = zeros(steps, N)
    separation_dist = zeros(steps)

    println("🔄 Ejecutando pasos temporales...")
    for k in 1:steps
        s_k = Float64(inputs[k])
        rz = 1.0 - 2.0 * s_k
        rx = 2.0 * sqrt(s_k * (1.0 - s_k))

        # 1. Inyección
        rho_A = inject_state_EraseWrite(rho_A, 0, rz, rx=rx)
        rho_B = inject_state_EraseWrite(rho_B, 0, rz, rx=rx)

        # 2. Evolución
        for _ in 1:n_substeps
            rho_A = step_rk4(rho_A, H_op, dt)
            rho_B = step_rk4(rho_B, H_op, dt)
        end


        #rho = apply_global_dephasing(rho, g, 'Z')


        # 3. Medida
        for i in 1:N
            p_obj = label_to_obj[z_labels[i]]
            hist_A[k, i] = real(get(rho_A, p_obj, 0.0))
            hist_B[k, i] = real(get(rho_B, p_obj, 0.0))
        end
        # --- OPCIONAL: Aplicar dephasing después de la medición ---
        # Esto modificará rho_A y rho_B para el paso k+1.
        #if measurement_strength > 0.0
        #    rho_A = apply_global_dephasing(rho_A, measurement_strength, "Z")
        #    rho_B = apply_global_dephasing(rho_B, measurement_strength, "Z")
        #end

        # 4. Cálculo de Separación
        if k > 1 && inputs[k] != inputs[k-1]
            dist = 0.0
            for i in 1:N; dist += (hist_A[k, i] - hist_A[k-1, i])^2; end
            separation_dist[k] = sqrt(dist)
        end
        if k % 10 == 0; print("\rProgreso: $k/$steps"); end
    end

    # ==========================================================================
    # 4. GUARDADO DE DATOS (JLD2)
    # ==========================================================================
    println("\n💾 Guardando datos en JLD2...")
    save(joinpath(Experiment_path, "validation_data.jld2"), Dict(
        "N" => N, "inputs" => inputs,
        "hist_A" => hist_A, "hist_B" => hist_B,
        "separation_dist" => separation_dist
    ))

    # ==========================================================================
    # 5. LLAMADA A LA FUNCIÓN DE VISUALIZACIÓN
    # ==========================================================================
    # Aquí llamamos a la función externa que limpia el código
    plot_and_save_validation(hist_A, hist_B, separation_dist, N, steps, Experiment_path)

    println("✅ ¡Proceso finalizado con éxito!")
end

run_nathan_validation_task()