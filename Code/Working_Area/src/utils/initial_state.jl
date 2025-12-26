# ==============================================================================
# src/utils/state_prep.jl
# Preparación de estados iniciales
# ==============================================================================
function initial_state_all_zeros(n_qubits::Int)
    rho = Operator()
    
    # Coeficiente c_k = Tr(rho P) / 2^N. 
    # Para <Z>=1, el coeficiente es 1/2^N.
    norm_coeff = (1.0 / 2.0)^n_qubits + 0.0im
    
    # --- CORRECCIÓN AQUÍ: PARÉNTESIS AÑADIDOS ---
    # Antes: 1 << n_qubits - 1  -> Julia hacía 1 << (6-1) = 32 (MAL)
    # Ahora: (1 << n_qubits) - 1 -> Julia hace 64 - 1 = 63 (BIEN)
    limit = (1 << n_qubits) - 1
    
    for z_mask_val in 0:limit
        p = PauliString(0, z_mask_val)
        rho[p] = norm_coeff
    end
    
    return rho
end





"""
    initial_state_x_first(n_qubits::Int)

Devuelve el estado ρ = X ⊗ I ⊗ ... ⊗ I.
Este es un estado "Sparse" (disperso) ideal para simulaciones grandes (N=20),
ya que solo contiene 1 término en lugar de 2^N.
"""
function initial_state_x_first(n_qubits::Int)
    rho = Operator()
    
    # Creamos una PauliString con X en el primer qubit (bit 0)
    # x_mask = 1 (binario ...001), z_mask = 0
    p_x1 = PauliString(1, 0) 
    
    # El coeficiente es 1.0
    rho[p_x1] = 1.0 + 0.0im
    
    return rho
end


