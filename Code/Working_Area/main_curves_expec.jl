using LinearAlgebra
using Statistics
using Plots
using Random 

# ==============================================================================
# 1. CARGA DE MÓDULOS
# ==============================================================================
include("src/operator_terms/pauli_algebra.jl")      
include("src/utils/dynamics.jl")
include("src/utils/injection_EraseWrite.jl")
include("src/utils/initial_state.jl") 
include("src/utils/quantum_channels.jl") 
include("src/utils/shot_noise.jl") 
include("src/utils/measurements.jl")
include("src/capacity_training/feature_extraction.jl") 
include("src/visualization/expectation_xj_vs_step.jl") 

# ==============================================================================
# 2. HAMILTONIANO TFIM (Correcto)
# ==============================================================================
function hamiltonian_tfim(n_qubits::Int, J_vec::Vector{Float64}, h::Float64)
    H_ops = Dict{PauliString, ComplexF64}()
    # -h * sum(X_i)
    for i in 0:(n_qubits-1)
        op_x = PauliString(1 << i, 0) 
        H_ops[op_x] = get(H_ops, op_x, 0.0im) - h
    end
    # +sum(J_i * Z_i * Z_{i+1})
    for i in 0:(n_qubits-2)
        mask_zz = (1 << i) | (1 << (i+1))
        op_zz = PauliString(0, mask_zz)
        H_ops[op_zz] = get(H_ops, op_zz, 0.0im) + J_vec[i+1]
    end
    return H_ops
end

# Función auxiliar de seguridad
function enforce_normalization!(rho::Operator, n_qubits::Int)
    # El coeficiente de la Identidad debe ser EXACTAMENTE 1/2^N
    # Si varía, todo explota.
    id_op = PauliString(0, 0)
    rho[id_op] = 1.0 / (2.0^n_qubits) + 0.0im
end

# ==============================================================================
# 3. REPLICACIÓN FIGURA 5 (ESTABLE)
# ==============================================================================
function run_figure_5_stable()
    println("🚀 REPLICANDO FIGURA 5 (Modo Exacto / Sin Truncamiento)...")

    # --- PARÁMETROS DEL PAPER ---
    N_QUBITS    = 6             
    NUM_STEPS   = 50            
    DT_TOTAL    = 10.0          # Paso gigante del paper
    N_SUBSTEPS  = 2000          # Integración muy fina
    
    J_SCALE     = 1.0           
    H_FIELD     = 10.0          
    G_STRENGTH  = 0.3           # Backaction fuerte
    
    NOISE_SHOTS = 1.5e6         

    # --- A. FÍSICA ---
    # Usamos J fijo o aleatorio (el paper usa aleatorio, pero probemos J=1 para ver difusión clara)
    # Jvec = (rand(N_QUBITS) .- 0.5) .* J_SCALE  <-- Lo del paper
    Jvec = ones(N_QUBITS) .* J_SCALE # <-- PRUEBA ESTO: J uniforme ayuda a ver la "ola" limpia
    
    println("⚙️ Hamiltoniano TFIM (h=$H_FIELD, J=1.0)...")
    H_evol = hamiltonian_tfim(N_QUBITS, Jvec, H_FIELD)
    
    rho = initial_state_all_zeros(N_QUBITS) 
    
    dt_rk4 = DT_TOTAL / N_SUBSTEPS
    scale_factor = 2.0^N_QUBITS
    
    # Base de Medición
    basis = build_reservoir_basis(N_QUBITS)
    offset_X = N_QUBITS # Los X están en la segunda mitad
    
    X_data = zeros(Float64, NUM_STEPS, N_QUBITS)

    println("⏳ Evolucionando (Sin Truncar)...")
    
    for k in 1:NUM_STEPS
        # --- B. INYECCIÓN ---
        # El paper rota de Z (0) a X (1).
        progress = k / NUM_STEPS
        r_x_val = sqrt(progress) 
        r_z_val = sqrt(max(0.0, 1.0 - r_x_val^2))
        
        rho = inject_state_EraseWrite(rho, 0, r_z_val, rx=r_x_val)

        # --- C. EVOLUCIÓN ESTABLE ---
        for _ in 1:N_SUBSTEPS
            rho = step_rk4(rho, H_evol, dt_rk4)
            # ¡IMPORTANTE! HEMOS QUITADO EL TRUNCATE.
            # truncate_operator!(rho, 2000) <- ESTO ERA EL CULPABLE
        end
        
        # Seguridad numérica
        enforce_normalization!(rho, N_QUBITS)

        # --- D. MEDICIÓN ---
        feats = extract_all_features(rho, basis)
        for q in 1:N_QUBITS
            raw_val = feats[offset_X + q]
            val_fisico = raw_val * scale_factor
            X_data[k, q] = apply_shot_noise(val_fisico, NOISE_SHOTS)
        end

        # --- E. DEPHASING ---
        g_effective = G_STRENGTH * DT_TOTAL
        rho = apply_global_dephasing(rho, g_effective, "X") # Base X
        
        if k % 5 == 0; print("\r   Progreso $k / $NUM_STEPS"); end
    end

    println("\n✅ Simulación completada.")
    plot_expectation_evolution(1:NUM_STEPS, X_data, N_QUBITS)
    println("💾 Gráfica guardada en Outputs/expectation_evolution.png")
end

run_figure_5_stable()