# src/operator_terms/hamiltonian.jl

using Random

struct PauliString
    x_mask::Int
    z_mask::Int
end

const Operator = Dict{PauliString, ComplexF64}



"""
    build_ising_hamiltonian(n_qubits::Int, J::Float64, h::Float64) -> Operator

Construye el Hamiltoniano del modelo de Ising de Campo Transverso (TFIM) **HOMOGÉNEO**.
Representa un sistema ideal ("cristal limpio") donde todas las interacciones son idénticas.

# Fórmula
H = -J ∑ Z_i Z_{i+1} - h ∑ X_i

# Parámetros
- `J`: Constante de acoplamiento. Es igual para todos los pares de vecinos.
- `h`: Campo magnético transverso. Es igual para todos los sitios.
"""

function build_ising_hamiltonian(n_qubits::Int, J::Float64, h::Float64)
    H = Operator()
    for i in 0:(n_qubits-2)
        z_mask = (1 << i) | (1 << (i+1))
        H[PauliString(0, z_mask)] = -J 
    end
    for i in 0:(n_qubits-1)
        H[PauliString(1 << i, 0)] = -h
    end
    return H
end




"""
    build_paper_hamiltonian(n_qubits, J_scale, h_const)
    Replica Eq (32): J aleatorio [-0.5, 0.5], h constante.
"""


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

function build_paper_hamiltonian(n_qubits::Int, J_scale::Float64, h::Float64; seed::Int=1234)
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
        H[p_x] = get(H, p_x, 0.0im) + h
    end
    
    return H
end