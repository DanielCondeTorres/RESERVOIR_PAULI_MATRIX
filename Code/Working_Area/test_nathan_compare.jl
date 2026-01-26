using LinearAlgebra
using Plots
using JLD2
using FileIO

# ==============================================================================
# 1. INCLUDES
# ==============================================================================
# Asegúrate de ejecutar este script desde la carpeta "Code/Working_Area"
# para que encuentre la carpeta "src" correctamente.
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
# 2. CONFIGURACIÓN ROBUSTA DE RUTAS Y PARÁMETROS
# ==============================================================================

# --- Rutas de Archivos (Corregido para evitar error "No file exists") ---
# Definimos la ruta base a "Input_Data" usando expanduser para que entienda el "~"
DIR_INPUT = expanduser("~/Desktop/IBM_EXAMENES/TFM/RESERVOIR_PAULI_MATRIX/Code/Input_Data")
NOMBRE_ARCHIVO = "6_3_2_all_zeros_12345.jld2"

# Construimos la ruta completa segura
INPUT_FILE = joinpath(DIR_INPUT, NOMBRE_ARCHIVO)

# Verificación de seguridad inmediata
if !isfile(INPUT_FILE)
    println("\n❌ ERROR CRÍTICO: No se encuentra el archivo de entrada.")
    println("   Buscando en: $INPUT_FILE")
    println("   Verifica que la carpeta 'Input_Data' esté en esa ruta.")
    exit()
end

# --- Parámetros de Simulación ---
N_SUBSTEPS = 100       

# Dónde inyectamos
BIT_PARA_CARGAR  = [1] 

# Qué medimos
INDICES_MEDICION = [2] 
PAULI_MATRIX     = "Z"    
DEPHASING_PAULI_MATRIX = "Z"

# Nombre del archivo de salida (se guardará donde ejecutes el script)
OUTPUT_FILE = "resultado_simulacion_$(PAULI_MATRIX)_$(INDICES_MEDICION).jld2"

# ==============================================================================
# 3. EJECUCIÓN
# ==============================================================================

function run_final_validation()
    println("🚀 Iniciando Simulación y Validación...")
    println("📂 Leyendo archivo: $NOMBRE_ARCHIVO")
    
    idx_injection = BIT_PARA_CARGAR[1]
    
    # --- A. Cargar Datos ---
    # 1. Usamos tu loader para sacar los params limpios (necesarios para la física)
    params = extract_metadata(INPUT_FILE)
    nathan_series = extract_nathan_data(INPUT_FILE, PAULI_MATRIX, INDICES_MEDICION)

    # 2. IMPORTANTE: Cargamos el JLD2 crudo para recuperar los metadatos ORIGINALES COMPLETOS
    # Esto asegura que al guardar, tengamos 'task', 'M', 'heisenberg', etc.
    datos_crudos = load(INPUT_FILE)
    meta_data_original_completa = datos_crudos["meta_data_dict"]
    
    println("\n⚠️ PARÁMETROS FÍSICOS:")
    println("   N_qubits = $(params["N_qubits"]) (Leído del loader)")
    println("   Metadatos originales recuperados: $(length(meta_data_original_completa)) claves.")

    # --- B. Preparar Simulación ---
    H_evol = hamiltonian_nathan_XX(params["N_qubits"], params["J"], params["h"])
    rho = initial_state_all_zeros(params["N_qubits"])
    
    observable_unico = create_pauli_observable(PAULI_MATRIX, INDICES_MEDICION)
    
    # LÓGICA DE ESCALA
    if INDICES_MEDICION == BIT_PARA_CARGAR
        factor_escala = 1.0
        println("   Escala: x1.0 (Inyección)")
    else
        factor_escala = 1.0 
        println("   Escala: x1.0 (Reservorio/Correlaciones)")
    end

    dt_small = params["dt"] / N_SUBSTEPS
    max_terms = 4^params["N_qubits"]
    
    my_history = zeros(Float64, params["num_steps"])

    # --- C. Bucle Temporal ---
    println("\n⚡ Calculando evolución temporal...")
    for k in 1:params["num_steps"]
        s_k = params["s_vec"][k]
        
        # 1. Inyección
        rz = 1.0 - 2.0 * s_k
        rx = 2.0 * sqrt(s_k * (1.0 - s_k))
        rho = inject_state_EraseWrite(rho, idx_injection, rz, rx=rx)
        
        # 2. Evolución
        for _ in 1:N_SUBSTEPS
            rho = step_rk4(rho, H_evol, dt_small)
            truncate_operator!(rho, max_terms)
        end
        
        # 3. Dephasing
        rho = apply_global_dephasing(rho, params["g"], DEPHASING_PAULI_MATRIX)
        
        # 4. Medición
        coeff, _ = measure_observable(rho, observable_unico, projective=false)
        my_history[k] = coeff * factor_escala 
        
        if k % 100 == 0
            print("\r   ⏳ Paso $k / $(params["num_steps"])")
        end
    end
    
    println("\n✅ Simulación Finalizada.")

    # ==========================================================================
    # --- D. GUARDADO EXACTO (REPLICANDO ESTRUCTURA) ---
    # ==========================================================================
    println("\n💾 Guardando archivo compatible...")

    # 1. Reconstruir la CLAVE COMPLEJA del diccionario (Vector de Tuplas)
    # Ejemplo: [(2, "Z")]
    clave_simulacion = []
    chars_pauli = split(PAULI_MATRIX, "")
    
    if length(chars_pauli) == length(INDICES_MEDICION)
        for (i, idx) in enumerate(INDICES_MEDICION)
            push!(clave_simulacion, (idx, string(chars_pauli[i])))
        end
    else
        for idx in INDICES_MEDICION
            push!(clave_simulacion, (idx, PAULI_MATRIX))
        end
    end
    
    # 2. Crear diccionario de resultados
    expect_dict_to_save = Dict()
    expect_dict_to_save[clave_simulacion] = my_history

    # 3. Guardar usando los metadatos ORIGINALES COMPLETOS
    save(OUTPUT_FILE, 
         "meta_data_dict", meta_data_original_completa, 
         "expect_dict", expect_dict_to_save)
         
    println("   ➜ Archivo guardado: $OUTPUT_FILE")
    println("   ➜ Validado: Contiene los metadatos completos originales.")

    # ==========================================================================

    # --- E. Graficar (Opcional) ---
    if !isnothing(nathan_series)
        len_comp = min(length(nathan_series), length(my_history), 200)
        err = (norm(nathan_series[1:len_comp] - my_history[1:len_comp]) / 
               norm(nathan_series[1:len_comp])) * 100
        
        println("\n📊 Generando gráfica comparativa (Error: $(round(err, digits=2))%)...")
        p = plot_validation_comparison(
            nathan_series, my_history, INDICES_MEDICION,     
            params["N_qubits"], err, PAULI_MATRIX, 200    
        )
        display(p)
    end
end

run_final_validation()