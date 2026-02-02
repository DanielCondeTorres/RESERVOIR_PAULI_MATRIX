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
h_val = 0.5
Experiment_name = "../Nathan_Foundation_AllToAll_FullZ_Operators"

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
include_rel("src/utils/injection_EraseWrite.jl")

# ==============================================================================
# 3. SIMULACIÓN BASADA EN OPERADORES (RK4 + COEFICIENTES)
# ==============================================================================
function run_nathan_task()
    println("Iniciando experimento (MODO OPERADORES): $Experiment_name")
    
    H_op = build_nathan_all_to_all_XX(N, h_val) 
    rho = initial_state_all_zeros(N) 
    
    Random.seed!(1234)
    inputs = rand(0:1, steps)

    println("Pre-calculando objetos PauliString...")
    
    # Necesitamos guardar z_labels por separado para la gráfica final
    z_labels = String[]
    label_to_obj = Dict{String, PauliString}()

    # 1. Individuales
    for i in 1:N
        label = join([k == i ? "Z" : "1" for k in 1:N])
        push!(z_labels, label)
        label_to_obj[label] = string_to_pauli_mask(label)
    end

    # 2. Correlaciones ZZ
    for i in 1:N, j in (i+1):N
        label = join([k == i || k == j ? "Z" : "1" for k in 1:N])
        label_to_obj[label] = string_to_pauli_mask(label)
    end

    all_labels = collect(keys(label_to_obj))
    expect_dict = Dict(lbl => zeros(Float64, steps) for lbl in all_labels)



    println(" Ejecutando pasos temporales con RK4...")
    for k in 1:steps
        s_k = Float64(inputs[k])

        # 1. INYECCIÓN
        rz = 1.0 - 2.0 * s_k
        rx = 2.0 * sqrt(s_k * (1.0 - s_k))
        rho = inject_state_EraseWrite(rho, 0, rz, rx=rx)
        
        
        #MEDICION PREVIA
        for label in all_labels
            p_obj = label_to_obj[label] 
            expect_dict[label][k] = real(get(rho, p_obj, 0.0))
        end
    # FIN MEDICION PREVIA
        # 2. EVOLUCIÓN
        for _ in 1:n_substeps
            rho = step_rk4(rho, H_op, dt)
        end
        # 3. MEDIDA CORREGIDA: Usamos el objeto ya traducido
        #for label in all_labels
        #    p_obj = label_to_obj[label] 
        #    expect_dict[label][k] = real(get(rho, p_obj, 0.0))
        #end
        if k % 10 == 0; print("\rStep $k/$steps"); end
    end

    # E. GUARDADO Y GRÁFICA
    ruta_experimento = save_qrc_results_jld2(N, steps, T_evol, h_val, inputs, expect_dict, Experiment_name)
    
    # Usamos las z_labels que guardamos arriba
    history_single_z = zeros(Float64, steps, N)
    for i in 1:N
        history_single_z[:, i] = expect_dict[z_labels[i]]
    end

    println("\n Generando gráfica...")
    plot_expectation_evolution_easy(1:steps, history_single_z, N, ruta_experimento, inputs)
    plot_zz_correlations(1:steps, expect_dict, ruta_experimento, inputs)
    println(" ¡Listo!")

        # 'history_single_z' es tu matriz de (steps, N)
    heatmap(1:steps, 1:N, history_single_z', 
    c=:viridis, 
    xlabel="Step (k)", 
    ylabel="Qubit Index",
    title="Reservoir Dynamics: Nathan Foundation",
    clabel="Expectation <Z>")
end

run_nathan_task()