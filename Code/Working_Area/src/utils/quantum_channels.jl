# src/utils/quantum_channels.jl
"""
    apply_global_dephasing(rho::Operator, g::Float64, axis::String)

Aplica el canal de desfase (medición débil) en la base indicada ('X', 'Y', o 'Z').
Implementa la backaction física:
- "Z": Penaliza términos que anticonmutan con Z (X, Y).
- "X": Penaliza términos que anticonmutan con X (Z, Y).
- "Y": Penaliza términos que anticonmutan con Y (X, Z).

Fórmula: coeff_new = coeff * exp(-g²/2)^n_anticonmuta
"""
function apply_global_dephasing(rho::Operator, g::Float64, axis::String)
    new_rho = Operator()
    decay_base = exp(-(g^2) / 2.0)
    axis_norm = uppercase(axis)
    
    # --- VARIABLES DE DEBUG ---
    total_terms = 0
    killed_terms = 0
    # --------------------------

    for (p, coeff) in rho
        if abs(coeff) < 1e-15; continue; end
        total_terms += 1
        
        n_anti = 0
        
        if axis_norm == "Z"
            # Si medimos Z, mueren los que tienen X (X e Y)
            n_anti = count_ones(p.x_mask)
            
        elseif axis_norm == "X"
            # Si medimos X, mueren los que tienen Z (Z e Y)
            n_anti = count_ones(p.z_mask)
            
        elseif axis_norm == "Y"
            n_anti = count_ones(p.x_mask ⊻ p.z_mask)
        end
        
        if n_anti > 0
            new_rho[p] = coeff * (decay_base ^ n_anti)
            killed_terms += 1
        else
            new_rho[p] = coeff
        end
    end

    # --- IMPRIMIR DEBUG ---
    # Esto imprimirá en cada paso qué está pasando.
    # Si ves que los números cambian al cambiar axis="X" a "Z", ¡funciona!
    if total_terms > 0
        println("🔎 DEBUG Dephasing [$axis_norm]: De $total_terms términos, se amortiguaron $killed_terms")
    end
    # ----------------------
    
    return new_rho
end


# Lógica por ejemplo para el dephasing en Y:
        # X (1,0) -> 1^0 = 1 (Anticonmuta con Y -> Muere)
        # Z (0,1) -> 0^1 = 1 (Anticonmuta con Y -> Muere)
        # Y (1,1) -> 1^1 = 0 (Conmuta con Y -> Vive)
        # I (0,0) -> 0^0 = 0 (Conmuta con Y -> Vive)




        # ==============================================================================
# DEPHASING PARA MATRICES DENSAS (Schrödinger)
# ==============================================================================


# --- FUNCIONES AUXILIARES DE ROTACIÓN ---

function rotate_density_matrix(rho, axis, n_qubits; inverse=false)
    axis_norm = uppercase(axis)
    if axis_norm == "Z"
        return rho # No hace falta rotar para Z
    end
    
    # Definir el gate de rotación local
    gate = if axis_norm == "X"
        # Hadamard mueve Z -> X
        1/sqrt(2) * [1.0 1.0; 1.0 -1.0]
    elseif axis_norm == "Y"
        # Rx(π/2) mueve Z -> Y
        # U = exp(-i π/4 X)
        1/sqrt(2) * [1.0 -1.0im; -1.0im 1.0]
    end
    
    if inverse; gate = gate'; end
    
    # Construir el operador global U = g ⊗ g ⊗ g...
    U_global = gate
    for i in 2:n_qubits
        U_global = kron(U_global, gate)
    end
    
    # Aplicar transformación adjunta: U * rho * U'
    return U_global * rho * U_global'
end





using LinearAlgebra

"""
    apply_global_dephasing_matrix(rho, g, axis)

Aplica dephasing global en cualquier eje (X, Y, Z) sobre una matriz densa.
Utiliza la lógica de decaimiento: coeff * exp(-g²/2)^n_diff
"""
function apply_global_dephasing_matrix(rho::Matrix{ComplexF64}, g::Float64, axis::String="Z")
    if g <= 1e-9; return rho; end
    
    dim = size(rho, 1)
    n_qubits = Int(log2(dim))
    decay_factor = exp(-(g^2) / 2.0) # Usamos la misma constante que en tu código de Paulis
    
    # 1. ROTAR AL EJE DESEADO
    # Para medir en X, rotamos la base con Hadamards. 
    # Para medir en Y, rotamos con Rx(π/2).
    rho_rotated = rotate_density_matrix(rho, axis, n_qubits, inverse=false)
    
    # 2. APLICAR DEPHASING (Lógica de Distancia de Hamming)
    # En la base rotada, el dephasing siempre se ve como un "Z-dephasing"
    new_rho = copy(rho_rotated)
    for c in 1:dim
        for r in 1:dim
            if r == c; continue; end
            
            n_diff = count_ones((r-1) ⊻ (c-1)) # Distancia de Hamming
            if n_diff > 0
                new_rho[r, c] *= (decay_factor ^ n_diff)
            end
        end
    end
    
    # 3. VOLVER A LA BASE ORIGINAL
    return rotate_density_matrix(new_rho, axis, n_qubits, inverse=true)
end