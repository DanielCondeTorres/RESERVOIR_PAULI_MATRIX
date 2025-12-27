"""
    build_reservoir_basis(n_qubits::Int)

Genera una base RICA de observables (Features) siguiendo el paper:
1. Locales: X, Y, Z
2. Correlaciones vecines: XX, YY, ZZ (y opcionalmente mixtas como XZ si quisieras)
"""
function build_reservoir_basis(n_qubits::Int)
    basis = PauliString[]
    
    # --- 1. LOCALES (Peso 1) ---
    for i in 0:(n_qubits-1)
        # Z (0, 1)
        push!(basis, PauliString(0, 1 << i))
        # X (1, 0)
        push!(basis, PauliString(1 << i, 0))
        # Y (1, 1) -> Mascara en X y en Z
        push!(basis, PauliString(1 << i, 1 << i))
    end
    
    # --- 2. CORRELACIONES VECINAS (Peso 2) ---
    for i in 0:(n_qubits-2)
        # Máscaras para los dos qubits (i e i+1)
        mask_pair = (1 << i) | (1 << (i+1))
        
        # ZZ: (0, mask)
        push!(basis, PauliString(0, mask_pair))
        
        # XX: (mask, 0)
        push!(basis, PauliString(mask_pair, 0))
        
        # YY: (mask, mask)
        push!(basis, PauliString(mask_pair, mask_pair))
    end
    
    return basis
end


function extract_all_features(rho::Operator, basis::Vector{PauliString})
    features = zeros(Float64, length(basis))
    for (i, P) in enumerate(basis)
        # Extraemos la parte real del valor esperado <P>
        features[i] = real(get(rho, P, 0.0im))
    end
    return features
end