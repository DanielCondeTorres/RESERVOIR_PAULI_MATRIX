# src/utils/quantum_channels.jl
function apply_global_dephasingZ(rho::Operator, g::Float64)
    new_rho = Operator()
    decay_base = exp(-(g^2) / 2.0)
    
    for (p, coeff) in rho
        n_non_commuting = count_ones(p.x_mask)  # Corregido: cuenta X e Y (anticommuta con Z)
        new_rho[p] = coeff * (decay_base ^ n_non_commuting)
    end
    return new_rho
end


function apply_global_dephasingX(rho::Operator, g::Float64)
    new_rho = Operator()
    decay_base = exp(-(g^2) / 2.0)
    
    for (p, coeff) in rho
        n_non_commuting = count_ones(p.z_mask)
        new_rho[p] = coeff * (decay_base ^ n_non_commuting)
    end
    return new_rho
end

"""
    apply_global_dephasingY(rho, g)
    Aplica ruido que PRESERVA la base Y.
    Penaliza (borra) los términos X y Z.
"""
function apply_global_dephasingY(rho::Operator, g::Float64)
    new_rho = Operator()
    decay_base = exp(-(g^2) / 2.0)
    
    for (p, coeff) in rho
        if abs(coeff) < 1e-15; continue; end
        
        # Lógica:
        # X (1,0) -> 1^0 = 1 (Anticonmuta con Y -> Muere)
        # Z (0,1) -> 0^1 = 1 (Anticonmuta con Y -> Muere)
        # Y (1,1) -> 1^1 = 0 (Conmuta con Y -> Vive)
        # I (0,0) -> 0^0 = 0 (Conmuta con Y -> Vive)
        n_kill = count_ones(p.x_mask ⊻ p.z_mask) # Operador XOR
        
        new_rho[p] = coeff * (decay_base ^ n_kill)
    end
    return new_rho
end