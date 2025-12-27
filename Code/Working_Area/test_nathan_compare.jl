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
include("src/visualization/plotting_comparission_vs_nathan.jl") 
# ==============================================================================
# 2. CONFIGURACIÓN DE ÍNDICES (CRÍTICO)
# ==============================================================================
INPUT_FILE = "6_3_2_all_zeros_12345.jld2" 
N_SUBSTEPS = 100       
# Queremos el Bit 1 (el segundo qubit):
BIT_PARA_CARGAR  = [1,3]  # Para extract_nathan_data -> Posición 2 ("1Z1...") 
PAULI_MATRIX = "Z"
DEPHASING_PAULI_MATRIX = "X"  # Para canal de dephasing global
NOISE_SHOTS = 1.5e6

# ==============================================================================
# 3. EJECUCIÓN
# ==============================================================================
BIT_PARA_SIMULAR = BIT_PARA_CARGAR .+ 1
function run_final_validation()
    println("🚀 Iniciando Validación Final...")

    # A. Cargar Metadatos y Datos de Referencia
    params = extract_metadata(INPUT_FILE)
    # Usamos BIT_PARA_CARGAR [1] para encontrar la clave "1Z1111"
    nathan_series = extract_nathan_data(INPUT_FILE, PAULI_MATRIX, BIT_PARA_CARGAR)
    
    if isnothing(nathan_series)
        error("❌ No se encontró el qubit en el archivo. Revisa BIT_PARA_CARGAR.")
    end

    # B. Preparar Simulación
    # NOTA: Si esto falla, intenta cambiar a hamiltonian_nathan_XX
    H_evol = hamiltonian_nathan_ZZ(params["N_qubits"], params["J"], params["h"])
    rho = initial_state_all_zeros(params["N_qubits"])
    
    # El observable para medir en la simulación (Bit 1)
    observable_unico = create_pauli_observable(PAULI_MATRIX, BIT_PARA_SIMULAR)
    
    scale_factor = 2.0^params["N_qubits"]
    dt_small = params["dt"] / N_SUBSTEPS
    my_history = zeros(Float64, params["num_steps"])

    # C. Bucle Temporal
    for k in 1:params["num_steps"]
        # Injection
        s_k = params["s_vec"][k]
        rz = 1.0 - 2.0 * s_k
        rx = 2.0 * sqrt(s_k * (1.0 - s_k))
        rho = inject_state_EraseWrite(rho, 0, rz, rx=rx) 
        # Evolution
        for _ in 1:N_SUBSTEPS
            rho = step_rk4(rho, H_evol, dt_small)
            truncate_operator!(rho, 2000)
        end
        
        # Dephasing
        rho = apply_global_dephasing(rho, params["g"],DEPHASING_PAULI_MATRIX) 
        # Medición: Obtener coeficiente y aplicar escala después del ruido
        # como en tu versión que "daba bien".
        coeff = real(get(rho, observable_unico, 0.0im))# REVISAR MEDICION



        my_history[k] = apply_shot_noise(coeff, NOISE_SHOTS) * scale_factor
        if k % 100 == 0; print("\r   ⏳ Paso $k / $(params["num_steps"]) ..."); end
    end
    # D. Graficar
    println("\n✅ Finalizado.")

    #limite_a_graficar = min(length(nathan_series), length(my_history))
    limite_a_graficar = 100 # si no pongo nada lo hace con el anterior
    err = (norm(nathan_series[1:limite_a_graficar] - my_history[1:limite_a_graficar]) / norm(nathan_series[1:limite_a_graficar])) * 100
    
    # 2. LLAMAR A LA FUNCIÓN DE VISUALIZACIÓN
    # Esta función se encarga de:
    #   - Crear la etiqueta correcta (ej: "IZZII") usando BIT_PARA_CARGAR
    #   - Poner los títulos y leyendas
    #   - Guardar automáticamente en la carpeta Outputs
    p = plot_validation_comparison(
        nathan_series,       # 1. Datos ref
        my_history,          # 2. Datos sim
        BIT_PARA_CARGAR,     # 3. Índices (ej: [1,2])
        params["N_qubits"],  # 4. Total qubits
        err,                 # 5. Error %
        PAULI_MATRIX,                 # 6. Carácter (Char)
        limite_a_graficar    # 7. Límite (Int)
    )
    display(p)
end

run_final_validation()