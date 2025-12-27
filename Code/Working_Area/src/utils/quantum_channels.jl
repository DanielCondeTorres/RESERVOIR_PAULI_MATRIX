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
    
    # Factor de decaimiento base (Eq. 15 del paper)
    decay_base = exp(-(g^2) / 2.0)
    
    # Normalizamos el input a mayúsculas para evitar errores
    axis_norm = uppercase(axis)
    
    # Pre-chequeo de seguridad
    if axis_norm ∉ ["X", "Y", "Z"]
        error("❌ Eje de dephasing desconocido: '$axis'. Usa 'X', 'Y' o 'Z'.")
    end

    for (p, coeff) in rho
        # Filtrado de términos nulos por eficiencia
        if abs(coeff) < 1e-15; continue; end
        
        n_anti = 0
        
        if axis_norm == "Z"
            # Anticonmuta con Z si tiene componente X (es decir, es X o Y)
            n_anti = count_ones(p.x_mask)
            
        elseif axis_norm == "X"
            # Anticonmuta con X si tiene componente Z (es decir, es Z o Y)
            n_anti = count_ones(p.z_mask)
            
        elseif axis_norm == "Y"
            # Anticonmuta con Y si es X o Z puros.
            # Y (1,1) conmuta. I (0,0) conmuta.
            # X (1,0) y Z (0,1) anticonmutan -> Usamos XOR (⊻)
            n_anti = count_ones(p.x_mask ⊻ p.z_mask)
        end
        
        # Aplicamos el decaimiento
        if n_anti > 0
            new_rho[p] = coeff * (decay_base ^ n_anti)
        else
            new_rho[p] = coeff
        end
    end
    
    return new_rho
end


# Lógica por ejemplo para el dephasing en Y:
        # X (1,0) -> 1^0 = 1 (Anticonmuta con Y -> Muere)
        # Z (0,1) -> 0^1 = 1 (Anticonmuta con Y -> Muere)
        # Y (1,1) -> 1^1 = 0 (Conmuta con Y -> Vive)
        # I (0,0) -> 0^0 = 0 (Conmuta con Y -> Vive)
