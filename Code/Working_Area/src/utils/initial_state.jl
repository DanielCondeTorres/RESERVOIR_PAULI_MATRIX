# ==============================================================================
# src/utils/state_prep.jl
# Preparación de estados iniciales
# ==============================================================================

"""
    initial_state_all_zeros(n_qubits::Int)

Devuelve el estado puro |000...0> (todos los qubits en |0>) en la base de Pauli.
Un término Pauli P tiene coeficiente (1/2)^N si P solo contiene I y Z.
En representación binaria, esto significa que x_mask debe ser 0.
"""
function initial_state_all_zeros(n_qubits::Int)
    rho = Operator()
    
    # Coeficiente de normalización para N qubits
    norm_coeff = (1.0 / 2.0)^n_qubits + 0.0im
    
    # Iteramos sobre todas las posibles combinaciones de I y Z
    # Esto corresponde a iterar z_mask desde 0 hasta 2^N - 1, manteniendo x_mask en 0.
    for z_mask_val in 0:(1 << n_qubits - 1)
        p = PauliString(0, z_mask_val) # x_mask es siempre 0
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