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