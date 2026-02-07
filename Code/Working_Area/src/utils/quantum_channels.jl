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
function apply_global_dephasing_schrodinger(rho::Matrix{ComplexF64}, g::Float64, axis::String="Z")
    # Si g es muy pequeño, no hacemos nada
    if g <= 1e-9; return rho; end
    
    # Obtenemos dimensión y número de qubits
    dim = size(rho, 1)
    # N no se pasa explícitamente, pero es log2(dim)
    
    # Factor de decaimiento por qubit anticonmutando.
    # Si usamos la lógica lineal: factor = exp(-g)
    decay_factor = exp(-g)
    
    # Copiamos para no mutar el original si no se desea
    new_rho = copy(rho)
    
    # IMPLEMENTACIÓN OPTIMIZADA POR DISTANCIA DE HAMMING
    # (Solo funciona nativamente para eje Z, que es lo que usas)
    if uppercase(axis) == "Z"
        for c in 1:dim
            for r in 1:dim
                if r == c; continue; end # La diagonal (poblaciones) no decae en dephasing Z
                
                # Truco binario:
                # Los índices en Julia son 1-based, restamos 1 para tener bits (0...2^N-1)
                val_r = r - 1
                val_c = c - 1
                
                # XOR (⊻) nos dice en qué bits (qubits) difieren el estado fila y col
                diff_bits = val_r ⊻ val_c
                
                # Contamos cuántos qubits son diferentes (Distancia Hamming)
                n_diff = count_ones(diff_bits)
                
                # Aplicamos el castigo: decae una vez por cada qubit diferente
                if n_diff > 0
                    new_rho[r, c] *= (decay_factor ^ n_diff)
                end
            end
        end
    else
        # Fallback simple si quisieras X o Y (aunque no lo usas ahora)
        # Sería rotar la base, aplicar Z, y desrotar.
        println("⚠️ Aviso: apply_global_dephasing denso optimizado solo para Z.")
    end
    
    return new_rho
end