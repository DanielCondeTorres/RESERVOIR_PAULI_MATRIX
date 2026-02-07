using LinearAlgebra
using Plots
using Random
using JLD2
using FileIO
using Dates
using Statistics
using StatsPlots
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
include_rel("src/utils/quantum_channels.jl")
include_rel("src/metrics/separability_metrics.jl")
include_rel("src/visualization/separability_plots.jl")
# ==============================================================================
# 3. FUNCIÓN PRINCIPAL DE SIMULACIÓN Y VALIDACIÓN
# ==============================================================================
function run_nathan_validation_task()
# ==========================================================================
# 1. PARÁMETROS "HIGH CONTRAST" (Optimizados para Separación)
# ==========================================================================
    N = 6
    steps = 100
    T_evol = 1.0 # Ajustado para mejor propagación
    h_val = 1.0 # OPTIMAL: Balance entre ESP y separabilidad
    n_substeps = 50 # Aumentado para precisión numérica
    dt = T_evol / n_substeps
    g = 0.05 # Aumentado ligeramente para fading sin matar señal
# Ruta
    Experiment_path = joinpath(SCRIPT_DIR, "../RESULTADOS_ESP_SEPARACION_FIX")
if !isdir(Experiment_path); mkpath(Experiment_path); end
println("🚀 Iniciando Simulación CORREGIDA")
# ==========================================================================
# 2. HAMILTONIANO DESORDENADO (Rompe la simetría)
# ==========================================================================
# Vector manual de 15 valores para N=6 (Mezcla fuertes y débiles)
    J_para_nathan = [1.23, 0.45, 1.89, 0.76, 1.12,
0.34, 1.56, 0.98, 1.41,
0.67, 1.21, 0.55,
1.78, 0.82,
1.05]
# ESCALADO: Dividimos por sqrt(N-1) para varianza finita como en modelos SK
    scale_factor = sqrt(N - 1)  # ≈2.236
    J_para_nathan = J_para_nathan / scale_factor
# Asegúrate de que tu función hamiltonian_nathan_XX acepte este vector
# Si no, usa build_random_all_to_all_XX
    H_op = hamiltonian_nathan_XX(N, J_para_nathan, h_val)
# ==========================================================================
# 3. PREPARACIÓN DE ESTRUCTURAS
# ==========================================================================
# Generamos TODAS las etiquetas (Z y ZZ)
    z_labels = [join([k == i ? "Z" : "1" for k in 1:N]) for i in 1:N]
    zz_labels = String[]
for i in 1:N, j in (i+1):N
        lbl = ["1" for _ in 1:N]; lbl[i] = "Z"; lbl[j] = "Z"
push!(zz_labels, join(lbl))
end
    all_labels = [z_labels; zz_labels]
    label_to_obj = Dict(lbl => string_to_pauli(lbl) for lbl in all_labels)
# Diccionarios para guardar la historia completa
    dict_A = Dict(lbl => zeros(steps) for lbl in all_labels)
    dict_B = Dict(lbl => zeros(steps) for lbl in all_labels)
# Matrices para compatibilidad con tus plots viejos
    hist_A = zeros(steps, N)
    hist_B = zeros(steps, N)
    separation_dist = zeros(steps)
    rho_A = initial_state_all_zeros(N)
    rho_B = initial_state_all_ones(N)
    Random.seed!(1234)
    inputs = rand(0:1, steps)
# ==========================================================================
# 4. BUCLE TEMPORAL (Aquí estaba el error)
# ==========================================================================
println("🔄 Ejecutando pasos temporales...")
for k in 1:steps
        s_k = Float64(inputs[k])
# --- CORRECCIÓN DE INYECCIÓN ---
# Inyectamos en Eje Y Puro (+Y o -Y) para mejor NS-ESP y separabilidad.
# Asumimos que inject_state_EraseWrite soporta ry (componente imaginaria para Y).
# Si no, modifica la función para psi = [rz_val, rx_val + 1im * ry_val]
# Para input 0: ry = +sin(pi/4), input 1: ry = -sin(pi/4)
# theta fijo en pi/4 para amplitud, pero signo en ry.
        theta = pi/4
        rz_val = cos(theta)  # √2/2 para ambos
        rx_val = 0.0  # Cambiado a 0 para eje Y
        ry_val = (s_k == 0 ? 1.0 : -1.0) * sin(theta)  # ± √2/2
# 1. Inyección (Usamos ry_val; ajusta la llamada si es necesario)
        rho_A = inject_state_EraseWrite(rho_A, 0, rz_val, rx=rx_val, ry=ry_val)  # Añade ry si la función lo permite
        rho_B = inject_state_EraseWrite(rho_B, 0, rz_val, rx=rx_val, ry=ry_val)
# 2. Evolución RK4
for _ in 1:n_substeps
            rho_A = step_rk4(rho_A, H_op, dt)
            rho_B = step_rk4(rho_B, H_op, dt)
end
# 3. Dephasing (DESCOMENTADO E IMPRESCINDIBLE)
        rho_A = apply_global_dephasing(rho_A, g, "Z")
        rho_B = apply_global_dephasing(rho_B, g, "Z")
# 4. Medida y LLENADO DE DICCIONARIOS (CRÍTICO)
for lbl in all_labels
            p_obj = label_to_obj[lbl]
# Medimos
            val_A = real(get(rho_A, p_obj, 0.0))
            val_B = real(get(rho_B, p_obj, 0.0))
# Guardamos en Diccionario (para los Boxplots y PCA)
            dict_A[lbl][k] = val_A
            dict_B[lbl][k] = val_B
# Guardamos en Matriz (para los plots viejos)
if count(c -> c == 'Z', lbl) == 1
                idx = findfirst('Z', lbl)
                hist_A[k, idx] = val_A
                hist_B[k, idx] = val_B
end
end
# 5. Cálculo de Separación Euclídea COMPLETA (Z + ZZ)
if k > 1 && inputs[k] != inputs[k-1]
            dist_sq = 0.0
for lbl in all_labels
# Sumamos diferencias al cuadrado de todos los observables
                dist_sq += (dict_A[lbl][k] - dict_A[lbl][k-1])^2
end
            separation_dist[k] = sqrt(dist_sq)
end
if k % 10 == 0; print("\rProgreso: $k/$steps"); end
end
# ==========================================================================
# 5. GENERACIÓN DE GRÁFICAS
# ==========================================================================
println("\n📊 Generando gráficas con datos llenos...")
# 1. Plots estilo Nathan (Validación ESP)
plot_and_save_validation_full(dict_A, dict_B, separation_dist, N, steps, Experiment_path)
# 2. Boxplots de Separación (Ahora sí saldrán datos)
plot_separability_boxplots(dict_A, inputs, N, Experiment_path)
# 3. Guardado
save(joinpath(Experiment_path, "validation_data_FIXED.jld2"), Dict(
"N" => N, "inputs" => inputs,
"dict_A" => dict_A, "separation" => separation_dist
    ))
println("✅ ¡LISTO! Revisa la carpeta RESULTADOS_ESP_SEPARACION_FIX")
end
# Ejecutar y valores bajos en la echo property...
run_nathan_validation_task()