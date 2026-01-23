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
# 2. CONFIGURACIÓN (¡CÁMBIALO A TU GUSTO!)
# ==============================================================================
INPUT_FILE = "6_3_2_all_zeros_12345.jld2" 
N_SUBSTEPS = 100       

# --- Dónde inyectamos ---
BIT_PARA_CARGAR  = [1] 

# --- Qué medimos ---
# PRUEBA AQUÍ LO QUE QUIERAS:
# - Qubit 1 solo: INDICES_MEDICION = [1], PAULI_MATRIX = "Z"
# - Qubit 2 solo: INDICES_MEDICION = [2], PAULI_MATRIX = "Z"
# - Correlación:  INDICES_MEDICION = [1, 2], PAULI_MATRIX = "ZZ"
INDICES_MEDICION = [2] 
PAULI_MATRIX     = "Z"    

DEPHASING_PAULI_MATRIX = "Z"

# ==============================================================================
# 3. EJECUCIÓN
# ==============================================================================

function run_final_validation()
    println("🚀 Iniciando Validación Inteligente...")
    
    idx_injection = BIT_PARA_CARGAR[1]
    
    # --- A. Cargar Datos ---
    params = extract_metadata(INPUT_FILE)
    nathan_series = extract_nathan_data(INPUT_FILE, PAULI_MATRIX, INDICES_MEDICION)
    
    println("\n⚠️ PARÁMETROS:")
    println("   N_qubits = $(params["N_qubits"])")
    println("   Inyección: Qubit $idx_injection")
    println("   Medición:  Qubits $INDICES_MEDICION ($PAULI_MATRIX)")

    # --- B. Preparar Simulación ---
    H_evol = hamiltonian_nathan_XX(params["N_qubits"], params["J"], params["h"])
    rho = initial_state_all_zeros(params["N_qubits"])
    
    observable_unico = create_pauli_observable(PAULI_MATRIX, INDICES_MEDICION)
    
    # ==========================================================================
    #  LÓGICA AUTOMÁTICA DE ESCALA
    # ==========================================================================
    # Si medimos EXACTAMENTE el qubit donde inyectamos, no escalamos (ya vale 1.0).
    # Si medimos cualquier otra cosa (vecinos o correlaciones), escalamos por 2^N.
    if INDICES_MEDICION == BIT_PARA_CARGAR
        factor_escala = 1.0
        println("Midiendo Inyección: Escala DESACTIVADA (x1.0)")
    else
        factor_escala = 1.0 #2.0^params["N_qubits"]
        println(" Midiendo Reservorio: Escala ACTIVADA (x$factor_escala)")
    end
    # ==========================================================================

    dt_small = params["dt"] / N_SUBSTEPS # DIVIDIR??
    max_terms = 4^params["N_qubits"]
    println("   Truncado: $max_terms (Full)")

    my_history = zeros(Float64, params["num_steps"])

    # --- C. Bucle Temporal ---
    for k in 1:params["num_steps"]
        s_k = params["s_vec"][k]
        
        # 1. INYECCIÓN
        rz = 1.0 - 2.0 * s_k
        rx = 2.0 * sqrt(s_k * (1.0 - s_k))
        rho = inject_state_EraseWrite(rho, idx_injection, rz, rx=rx)
        
        # 2. EVOLUCIÓN
        for _ in 1:N_SUBSTEPS
            rho = step_rk4(rho, H_evol, dt_small)
            truncate_operator!(rho, max_terms) # ¡IMPORTANTE: 4100!
        end
        
        # 3. DEPHASING
        rho = apply_global_dephasing(rho, params["g"], DEPHASING_PAULI_MATRIX)
        
        # 4. MEDICIÓN CON ESCALA AUTOMÁTICA
        #coeff = real(get(rho, observable_unico, 0.0im))
        coeff, _ = measure_observable(rho, observable_unico, projective=false)
        # Aquí aplicamos el factor que calculamos arriba
        my_history[k] = coeff * factor_escala 
        
        if k % 100 == 0
            print("\r   ⏳ Paso $k / $(params["num_steps"]) ...")
        end
    end
    
    println("\n\nSimulación Finalizada.")

    # --- D. Graficar ---
    if !isnothing(nathan_series)
        limite_a_graficar = 200
        len_comp = min(length(nathan_series), length(my_history), limite_a_graficar)
        err = (norm(nathan_series[1:len_comp] - my_history[1:len_comp]) / 
               norm(nathan_series[1:len_comp])) * 100
        
        println("   Error relativo: $(round(err, digits=2))%")
        
        p = plot_validation_comparison(
            nathan_series, my_history, INDICES_MEDICION,     
            params["N_qubits"], err, PAULI_MATRIX, limite_a_graficar    
        )
        display(p)
    else
        p = plot(my_history[1:200], label="Simulación", lw=2, color=:blue,
                 title="Medición $PAULI_MATRIX en $INDICES_MEDICION", xlabel="Pasos", ylabel="Valor")
        display(p)
    end
end

run_final_validation()