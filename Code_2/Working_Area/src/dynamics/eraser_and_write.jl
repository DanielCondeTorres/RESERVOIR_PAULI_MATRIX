function inject_state_EraseWrite_pauli(rho::Operator, qubit_idx::Int, rz::Float64, rx::Float64=0.0, ry::Float64=0.0)
    new_rho = Operator()
    bit_loc = 1 << qubit_idx
    
    for (p, coeff) in rho
        # ERASE: Solo sobreviven términos con Identidad en qubit_idx
        is_identity_at_k = ((p.x_mask >> qubit_idx) & 1 == 0) && ((p.z_mask >> qubit_idx) & 1 == 0)
        
        if is_identity_at_k
            # WRITE: El término Identidad se mantiene, y se expanden X, Y, Z
            new_rho[p] = get(new_rho, p, 0.0im) + coeff
            if abs(rz) > 1e-15
                pz = PauliString(p.x_mask, p.z_mask | bit_loc)
                new_rho[pz] = get(new_rho, pz, 0.0im) + coeff * rz
            end
            if abs(rx) > 1e-15
                px = PauliString(p.x_mask | bit_loc, p.z_mask)
                new_rho[px] = get(new_rho, px, 0.0im) + coeff * rx
            end
            if abs(ry) > 1e-15
                py = PauliString(p.x_mask | bit_loc, p.z_mask | bit_loc)
                new_rho[py] = get(new_rho, py, 0.0im) + coeff * ry
            end
        end
    end
    return new_rho
end