function include_rel(path...)
    include(joinpath(SCRIPT_DIR, path...))
end
using Random
struct PauliString
    x_mask::Int
    z_mask::Int
end

"""
    label_to_pauli(lbl::String, n_qubits::Int) -> PauliString

Convierte una etiqueta de texto (ej: "Z11", "1X1", "YY1") en un objeto PauliString.
- '1' representa la Identidad (no activa ningún bit).
- 'X' activa el bit correspondiente en x_mask.
- 'Z' activa el bit correspondiente en z_mask.
- 'Y' activa los bits en AMBAS máscaras (X e Z).
"""
function label_to_pauli(lbl::String, n_qubits::Int)
    x_mask = 0
    z_mask = 0
    
    # Recorremos cada carácter de la etiqueta
    for (i, char) in enumerate(lbl)
        # i-1 porque en Julia las máscaras de bits suelen empezar en 0
        if char == 'X'
            x_mask |= (1 << (i-1))
        elseif char == 'Z'
            z_mask |= (1 << (i-1))
        elseif char == 'Y'
            # Y se representa como X=1 y Z=1 en el mismo qubit
            x_mask |= (1 << (i-1))
            z_mask |= (1 << (i-1))
        end
        # Si es '1', no se hace nada (permanece como identidad local)
    end
    
    return PauliString(x_mask, z_mask)
end
function generate_obs_labels(n_qubits::Int)
    labels = String[]
    # Iteramos por las tres bases principales
    for b_char in ["X", "Y", "Z"]
        # 1-body: Un operador en un qubit, el resto identidad (1)
        for i in 1:n_qubits
            lbl = ["1" for _ in 1:n_qubits]
            lbl[i] = b_char
            push!(labels, join(lbl))
        end
        # 2-body: Dos operadores iguales, el resto identidad
        for i in 1:n_qubits, j in (i+1):n_qubits
            lbl = ["1" for _ in 1:n_qubits]
            lbl[i] = b_char; lbl[j] = b_char
            push!(labels, join(lbl))
        end
        # 3-body (Tríadas): Tres operadores iguales, el resto identidad
        for i in 1:n_qubits, j in (i+1):n_qubits, k in (j+1):n_qubits
            lbl = ["1" for _ in 1:n_qubits]
            lbl[i] = b_char; lbl[j] = b_char; lbl[k] = b_char
            push!(labels, join(lbl))
        end
    end
    return labels
end
const Operator = Dict{PauliString, ComplexF64}
# ==============================================================================
# 1. ESTRUCTURAS Y ÁLGEBRA BÁSICA
# ==============================================================================
struct PauliString
    x_mask::Int
    z_mask::Int
end
const Operator = Dict{PauliString, ComplexF64}

# Multiplicación: P1 * P2 = (fase) * P3
function multiply_paulis(p1::PauliString, p2::PauliString)
    new_x = p1.x_mask ⊻ p2.x_mask
    new_z = p1.z_mask ⊻ p2.z_mask
    phase = 1.0 + 0.0im
    combined = p1.x_mask | p1.z_mask | p2.x_mask | p2.z_mask
    n_bits = combined == 0 ? 0 : floor(Int, log2(combined)) + 1
    for i in 0:n_bits
        x1, z1 = (p1.x_mask >> i) & 1, (p1.z_mask >> i) & 1
        x2, z2 = (p2.x_mask >> i) & 1, (p2.z_mask >> i) & 1
        t1, t2 = x1 + 2*z1, x2 + 2*z2
        if t1 == 0 || t2 == 0 || t1 == t2; continue; end
        if (t1==1 && t2==2); phase *= im;      # XY = iZ
        elseif (t1==2 && t2==1); phase *= -im; 
        elseif (t1==2 && t2==3); phase *= im;  # YZ = iX
        elseif (t1==3 && t2==2); phase *= -im;
        elseif (t1==3 && t2==1); phase *= im;  # ZX = iY
        elseif (t1==1 && t2==3); phase *= -im;
        end
    end
    return PauliString(new_x, new_z), phase
end

# --- FUNCIÓN DE APOYO: Suma de Operadores ---
function add_ops(O::Operator, K::Operator, factor::ComplexF64)
    res = copy(O)
    for (p, coeff) in K
        res[p] = get(res, p, 0.0im) + coeff * factor
    end
    return res
end
