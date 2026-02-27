function build_nathan_all_to_all_XX(n_qubits::Int, h::Float64, J_vec::Union{Vector{Float64}, Nothing}=nothing)
    n_pairs = Int(n_qubits * (n_qubits - 1) / 2)
    H = Operator()
    
    # --- LÓGICA DE SELECCIÓN DE J ---
    # Definimos Js explícitamente para claridad
    Js = 1.0 
    
    actual_J = if isnothing(J_vec) || length(J_vec) != n_pairs
        if isnothing(J_vec)
            println("🎲 J_vec is 'nothing'. Generating random J_ij en [-1, 1] (Js=$Js) para N=$n_qubits...")
        else
            println("⚠️ WARNING: J_vec length mismatch. Expected $n_pairs, got $(length(J_vec)).")
            println("🎲 Generating new random J_ij values...")
        end
        # Genera valores uniformes en [-Js, Js]
        # (rand(n_pairs) .- 0.5) * 2.0 * Js  <-- Esto da el rango [-1, 1] si Js=1
        (rand(n_pairs) .- 0.5)*Js
    else
        J_vec
    end

    # 1. Campo Transversal: -(h/2) * Z_j
    # Según Mujal, h=10Js y se aplica como h/2
    for j in 0:(n_qubits - 1)
        p_z = PauliString(0, 1 << j) 
        H[p_z] = get(H, p_z, 0.0im) + h / 2.0
    end

    # 2. Interacción All-to-All: J_ij * X_i * X_j
    idx = 1 
    for i in 0:(n_qubits - 1)
        for j in (i + 1):(n_qubits - 1)
            J_ij = actual_J[idx]
            p_xx = PauliString((1 << i) | (1 << j), 0)
            H[p_xx] = get(H, p_xx, 0.0im) + J_ij
            idx += 1
        end
    end
    return H
end