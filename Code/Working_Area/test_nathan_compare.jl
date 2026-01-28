using LinearAlgebra
using Plots
using JLD2
using FileIO

# ==============================================================================
# 1. INCLUDES
# ==============================================================================
include("src/operator_terms/pauli_algebra.jl")      
include("src/utils/dynamics.jl")
include("src/utils/injection_EraseWrite.jl")
include("src/utils/initial_state.jl") 
include("src/utils/quantum_channels.jl") 
include("src/utils/shot_noise.jl") 
include("src/operator_terms/hamiltonian.jl") 
include("src/utils/measurements.jl")
include("src/loaders/data_loader.jl") 
include("src/visualization/plotting_time_series_vs_expectation_value.jl") 

# ==============================================================================
# 2. CONFIGURACIÓN
# ==============================================================================
DIR_DATOS = expanduser("~/Desktop/IBM_EXAMENES/TFM/RESERVOIR_PAULI_MATRIX/Code/Input_Data")
NOMBRE_ARCHIVO = "6_3_2_all_zeros_12345.jld2"
INPUT_FILE = joinpath(DIR_DATOS, NOMBRE_ARCHIVO)

if !isfile(INPUT_FILE)
    println("❌ Error: No se encuentra $INPUT_FILE")
    exit()
end

# --- CONTROL FÍSICO ---
ENABLE_INJECTION = true
ENABLE_DEPHASING = true

# Parámetros
N_SUBSTEPS = 100       
BIT_PARA_CARGAR  = [1] 
PAULI_MATRIX     = "Z"    
DEPHASING_PAULI_MATRIX = "Z"

# Ruta de salida
OUTPUT_FILE = expanduser("~/Desktop/IBM_EXAMENES/TFM/RESERVOIR_PAULI_MATRIX/Code/To_Nahan_ZZ/resultado_CON_ZZ_$(PAULI_MATRIX).jld2")
mkpath(dirname(OUTPUT_FILE))

# ==============================================================================
# 3. EJECUCIÓN
# ==============================================================================

function run_final_validation()
    println("🚀 Simulación con Correlaciones ZZ (Singles + Pairs)")
    
    idx_injection = BIT_PARA_CARGAR[1]
    
    # --- A. CARGA ---
    params = extract_metadata(INPUT_FILE)
    diccionario_completo_original = load(INPUT_FILE)
    meta_data_para_guardar = diccionario_completo_original["meta_data_dict"]
    
    N = params["N_qubits"]
    println("   N_qubits: $N")
    
    # --- B. PREPARACIÓN ---
    H_evol = hamiltonian_nathan_XX(N, params["J"], params["h"])
    rho = initial_state_all_zeros(N)
    
    # ---------------------------------------------------------
    # B.1 OBSERVABLES INDIVIDUALES (Z)
    # ---------------------------------------------------------
    lista_singles = []
    for i in 1:N
        push!(lista_singles, create_pauli_observable("Z", [i]))
    end
    
    # ---------------------------------------------------------
    # B.2 OBSERVABLES DE PARES (ZZ)
    # ---------------------------------------------------------
    # Generamos todas las combinaciones únicas (i, j) donde i < j
    lista_pares_ops = []
    lista_pares_idx = [] # Guardamos los índices (i,j) para luego saber quiénes son
    
    for i in 1:N
        for j in (i+1):N
            # Crea operador ZZ en posiciones i, j
            op_zz = create_pauli_observable("ZZ", [i, j])
            push!(lista_pares_ops, op_zz)
            push!(lista_pares_idx, (i, j))
        end
    end
    
    num_pares = length(lista_pares_ops)
    println("   Observables a medir: $N singles + $num_pares pares ZZ")

    # Matrices para guardar historia
    historia_singles = zeros(Float64, params["num_steps"], N)
    historia_pares   = zeros(Float64, params["num_steps"], num_pares)
    
    dt_small = params["dt"] / N_SUBSTEPS
    max_terms = 4^N

    # --- C. BUCLE TEMPORAL ---
    println("⚡ Calculando...")
    for k in 1:params["num_steps"]
        s_k = params["s_vec"][k]
        
        # 1. Inyección
        if ENABLE_INJECTION
            rz = 1.0 - 2.0 * s_k
            rx = 2.0 * sqrt(s_k * (1.0 - s_k))
            rho = inject_state_EraseWrite(rho, idx_injection, rz, rx=rx)
        end
        
        # 2. Evolución
        for _ in 1:N_SUBSTEPS
            rho = step_rk4(rho, H_evol, dt_small)
            truncate_operator!(rho, max_terms)
        end
        
        # 3. Dephasing
        if ENABLE_DEPHASING
            rho = apply_global_dephasing(rho, params["g"], DEPHASING_PAULI_MATRIX)
        end
        
        # 4. MEDICIONES
        
        # a) Medir Singles (Z)
        for i in 1:N
            val, _ = measure_observable(rho, lista_singles[i], projective=false)
            historia_singles[k, i] = val
        end
        
        # b) Medir Pares (ZZ)
        for p in 1:num_pares
            val_zz, _ = measure_observable(rho, lista_pares_ops[p], projective=false)
            historia_pares[k, p] = val_zz
        end
        
        if k % 100 == 0; print("\r   ⏳ Paso $k..."); end
    end
    println("\n✅ Terminado.")

    # ==========================================================================
    # --- D. GUARDADO (FORMATO STRING: "Z11...", "ZZ1...") ---
    # ==========================================================================
    println("💾 Guardando diccionario extendido...")
    expect_dict_nuevo = Dict{String, Any}()
    
    # 1. GUARDAR SINGLES (Z)
    for i in 1:N
        # Construir string "11Z111"
        chars_clave = collect(repeat("1", N))
        chars_clave[i] = 'Z'
        clave_str = join(chars_clave)
        
        expect_dict_nuevo[clave_str] = historia_singles[:, i]
    end
    
    # 2. GUARDAR PARES (ZZ)
    for p in 1:num_pares
        (i, j) = lista_pares_idx[p]
        
        # Construir string "1Z1Z11" (ejemplo i=2, j=4)
        chars_clave = collect(repeat("1", N))
        chars_clave[i] = 'Z'
        chars_clave[j] = 'Z'
        clave_str = join(chars_clave)
        
        expect_dict_nuevo[clave_str] = historia_pares[:, p]
        
        # Debug para ver si las claves salen bien
        if p <= 3; println("   Clave par generada: $clave_str"); end
    end

    save(OUTPUT_FILE, 
         "meta_data_dict", meta_data_para_guardar, 
         "expect_dict", expect_dict_nuevo)
         
    println("   ➜ Archivo: $OUTPUT_FILE")
    println("   ➜ Total claves guardadas: $(length(expect_dict_nuevo))")
    
    # Plot rápido de comparación (Single vs Correlation)
    # Graficamos el qubit 1 y la correlación entre 1 y 2
    if N >= 2
        p = plot(historia_singles[:, 1], label="Single <Z1>", lw=2)
        # Buscamos el índice del par (1,2)
        idx_pair_12 = findfirst(x -> x == (1, 2), lista_pares_idx)
        if !isnothing(idx_pair_12)
            plot!(p, historia_pares[:, idx_pair_12], label="Corr <Z1 Z2>", lw=2, linestyle=:dash)
        end
        display(p)
    end
end

run_final_validation()