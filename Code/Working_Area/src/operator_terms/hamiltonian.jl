# src/operator_terms/hamiltonian.jl

using Random

struct PauliString
    x_mask::Int
    z_mask::Int
end

const Operator = Dict{PauliString, ComplexF64}


"""
    build_paper_hamiltonian(n_qubits, J_scale, h_const) -> Operator

Construye un Hamiltoniano de Ising **DESORDENADO** (Eq. 32 del paper/réplica).
Representa un sistema con impurezas o tipo vidrio de espín (Spin Glass) en 1D.

# Fórmula
H = - ∑ J_i Z_i Z_{i+1} - h_{const} ∑ X_i

# Diferencia Clave
- `J_i`: Variable aleatoria uniforme centrada en 0. Rango: [-J_scale/2, +J_scale/2].
  Esto significa que algunos enlaces pueden ser ferromagnéticos (+) y otros antiferromagnéticos (-).
"""

function build_ising_hamiltonianZZ(n_qubits::Int, J_scale::Float64, h::Float64; seed::Int=1234)
    Random.seed!(seed) # 1. Reproducibilidad
    H = Operator()  
    # Término de interacción (J_ij aleatorio centrado en 0)
    for i in 0:(n_qubits-2)
        # 2. Distribución correcta según notas: [-J/2, J/2]
        J_val = J_scale * (rand() - 0.5)
        # 3. Construcción segura (+) y sin signo hardcodeado
        p_zz = PauliString(0, (1 << i) | (1 << (i+1)))
        H[p_zz] = get(H, p_zz, 0.0im) + J_val
    end
    # Término de campo transversal
    for i in 0:(n_qubits-1)
        p_x = PauliString(1 << i, 0)
        H[p_x] = get(H, p_x, 0.0im) - h
    end
    return H
end


"""
    build_paper_hamiltonian_nathan(n_qubits, J_scale, h) -> Operator

Replica EXACTAMENTE la Eq (32) del documento de Nathan Keenan:
H = - (h/2) ∑ Z_j + ∑ J_j X_j X_{j+1}
"""
function build_ising_hamiltonianXX(n_qubits::Int, J_scale::Float64, h::Float64; seed::Int=12345)
    Random.seed!(seed) 
    H = Operator()
    # 1. Término de Interacción: J_j * X_j * X_{j+1}
    # Nota: Nathan usa X para la interacción 
    for i in 0:(n_qubits-2)
        J_val = J_scale * (rand() - 0.5)
        # x_mask activo en los bits i e i+1
        p_xx = PauliString((1 << i) | (1 << (i+1)), 0)
        H[p_xx] = get(H, p_xx, 0.0im) + J_val
    end
    # 2. Término de Campo Transversal: -(h/2) * Z_j
    # Nota: Nathan usa Z para el campo y un factor 1/2 
    for i in 0:(n_qubits-1)
        p_z = PauliString(0, 1 << i)
        H[p_z] = get(H, p_z, 0.0im) - (h / 2.0)
    end
    return H
end










# AÑADIR ESTO AL FINAL DE src/operator_terms/hamiltonian.jl
# AÑADIR AL FINAL DE src/operator_terms/hamiltonian.jl

"""
    build_nathan_exact_hamiltonian(n_qubits::Int, J_vec::Vector{Float64}, h::Float64) -> Operator

Construye el Hamiltoniano EXACTO basado en la captura 10.47.09:
- Interacción: XX (Strings "X") entre vecinos.
- Signo J: POSITIVO (+J), como se ve en el código 'H += J...'.
- Campo: Z (Para permitir dinámica con interacción XX).
- J_vec: Vector exacto de acoplamientos cargado del archivo.
"""
function build_nathan_exact_hamiltonian(n_qubits::Int, J_vec::Vector{Float64}, h::Float64)
    H = Operator()
    # 1. Interacción XX con signo POSITIVO (+J)
    for i in 0:(n_qubits-2)
        # PauliString(x_mask, z_mask) -> XX tiene mascara X en i, i+1
        mask_xx = (1 << i) | (1 << (i+1))
        p_xx = PauliString(mask_xx, 0) 
        # ¡IMPORTANTE! Signo +J según la captura
        H[p_xx] = get(H, p_xx, 0.0im) + J_vec[i+1] 
    end
    # 2. Campo Transversal Z (Signo estándar -h/2)
    for i in 0:(n_qubits-1)
        p_z = PauliString(0, 1 << i) # Máscara Z
        H[p_z] = get(H, p_z, 0.0im) - (h / 2.0)
    end
    
    return H
end





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