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








"""
    build_nathan_all_to_all_XX(n_qubits::Int, h::Float64; seed::Int=12345) -> Operator

Construye el Hamiltoniano de la Eq (32) versión All-to-All:
H = - (h/2) ∑ Z_j + ∑_{i<j} J_{ij} X_i X_j

# Parámetros
- `n_qubits`: Número de sitios (ej. 6).
- `h`: Intensidad del campo magnético.
- `seed`: Para que los acoplamientos J_ij aleatorios sean siempre los mismos.
"""
function build_nathan_all_to_all_XX(n_qubits::Int, h::Float64; seed::Int=12345)
    Random.seed!(seed)
    H = Operator()

    # 1. Término de Campo Transversal: -(h/2) * Z_j en cada sitio
    # En tu máscara: PauliString(x_mask, z_mask)
    for j in 0:(n_qubits - 1)
        p_z = PauliString(0, 1 << j) # Z en el sitio j
        H[p_z] = get(H, p_z, 0.0im) - (h / 2.0)
    end

    # 2. Término de Interacción All-to-All: J_ij * X_i * X_j
    # Conectamos todos los pares posibles (i < j)
    for i in 0:(n_qubits - 1)
        for j in (i + 1):(n_qubits - 1)
            # J_ij aleatorio entre -0.5 y 0.5 (como indica Nathan en sus notas)
            J_ij = rand() - 0.5
            
            # X_i * X_j -> x_mask activo en bits i y j
            p_xx = PauliString((1 << i) | (1 << j), 0)
            
            # Sumamos al operador (por si acaso ya existiera el término)
            H[p_xx] = get(H, p_xx, 0.0im) + J_ij
        end
    end

    return H
end