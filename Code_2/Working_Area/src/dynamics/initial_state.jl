function initial_state_all_zeros(n_qubits::Int)
    rho = Operator()
    # Coeficiente c_k = Tr(rho P) / 2^N. 
    # Para <Z>=1, el coeficiente es 1/2^N.
    norm_coeff = (1.0 / 2.0)^n_qubits + 0.0im
    # --- CORRECCIÓN AQUÍ: PARÉNTESIS AÑADIDOS ---
    limit = (1 << n_qubits) - 1
    for z_mask_val in 0:limit
        p = PauliString(0, z_mask_val)
        rho[p] = norm_coeff
    end 
    return rho
end
function initial_state_all_ones(n_qubits::Int)
    rho = Operator()
    # El coeficiente base es 1/2^N
    base_coeff = (1.0 / 2.0)^n_qubits + 0.0im
    # Recorremos todas las combinaciones posibles de Z (2^N)
    limit = (1 << n_qubits) - 1
    for z_mask in 0:limit
        # Contamos cuántos Z hay en esta máscara (cuántos bits a 1)
        num_zs = count_ones(z_mask)
        sign = (num_zs % 2 == 0) ? 1.0 : -1.0
        p = PauliString(0, z_mask)
        rho[p] = sign * base_coeff
    end
    return rho
end