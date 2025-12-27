# src/operator_terms/hamiltonian.jl

using Random

struct PauliString
    x_mask::Int
    z_mask::Int
end

const Operator = Dict{PauliString, ComplexF64}


#######FUNCIONES HAMILTONIANO DE Nathan

"""
    hamiltonian_nathan_ZZ(n_qubits, J_exactos, h_exacto)

Construye el modelo Ising Estándar (Interacción ZZ, Campo X).
Corresponde a tu primer snippet de código.

Fórmula: H = - ∑ J_i Z_i Z_{i+1} - h ∑ X_i
"""
function hamiltonian_nathan_ZZ(n_qubits::Int, J_exactos::Vector{Float64}, h_exacto::Float64)
    H = Operator()
    
    # 1. Interacción entre vecinos (Z - Z)
    # Eje Z -> Segundo argumento de PauliString
    # Signo -> Negativo (-J)
    for i in 0:(n_qubits-2)
        z_mask = (1 << i) | (1 << (i+1))
        p_zz = PauliString(0, z_mask)
        H[p_zz] = -J_exactos[i+1]
    end

    # 2. Campo Transverso (X)
    # Eje X -> Primer argumento de PauliString
    # Signo -> Negativo (-h)
    for i in 0:(n_qubits-1)
        x_mask = (1 << i)
        p_x = PauliString(x_mask, 0)
        H[p_x] = -h_exacto
    end
    
    return H
end

"""
    hamiltonian_nathan_XX(n_qubits, J_exactos, h_exacto)

Construye el modelo Rotado/Dual (Interacción XX, Campo Z).
Corresponde a la captura del paper de Nathan (Eq. 32).

Fórmula: H = + ∑ J_i X_i X_{i+1} - (h/2) ∑ Z_i
"""
function hamiltonian_nathan_XX(n_qubits::Int, J_exactos::Vector{Float64}, h_exacto::Float64)
    H = Operator()
    
    # 1. Interacción entre vecinos (X - X)
    # Eje X -> Primer argumento de PauliString
    # Signo -> POSITIVO (+J), según la captura de Nathan
    for i in 0:(n_qubits-2)
        x_mask = (1 << i) | (1 << (i+1))
        p_xx = PauliString(x_mask, 0)
        H[p_xx] = +J_exactos[i+1] 
    end

    # 2. Campo Transverso (Z)
    # Eje Z -> Segundo argumento de PauliString
    # Signo -> Negativo pero a la MITAD (-h/2)
    for i in 0:(n_qubits-1)
        z_mask = (1 << i)
        p_z = PauliString(0, z_mask)
        H[p_z] = -(h_exacto / 2.0)
    end
    
    return H
end