# ==========================================
# ESTRUCTURAS Y ÁLGEBRA DE PAULI
# ==========================================

struct PauliString
    x_mask::Int
    z_mask::Int
end

# Alias para el Operador (Diccionario)
const Operator = Dict{PauliString, ComplexF64}

# Visualización bonita
function Base.show(io::IO, p::PauliString)
    print(io, "P($(string(p.x_mask, base=2)), $(string(p.z_mask, base=2)))")
end

# Multiplicación exacta bit a bit (Garantiza fase correcta)
function multiply_paulis(p1::PauliString, p2::PauliString)
    new_x = p1.x_mask ⊻ p2.x_mask
    new_z = p1.z_mask ⊻ p2.z_mask
    phase = 1.0 + 0.0im
    
    active_bits = (p1.x_mask | p1.z_mask) & (p2.x_mask | p2.z_mask)
    idx = 0
    while active_bits > 0
        if (active_bits & 1) == 1
            type1 = ((p1.x_mask >> idx) & 1) + 2 * ((p1.z_mask >> idx) & 1)
            type2 = ((p2.x_mask >> idx) & 1) + 2 * ((p2.z_mask >> idx) & 1)
            
            if type1 == 1 && type2 == 2     # XZ = -iY
                phase *= -1.0im
            elseif type1 == 2 && type2 == 1 # ZX = iY
                phase *= 1.0im
            elseif type1 == 2 && type2 == 3 # ZY = -iX
                phase *= -1.0im
            elseif type1 == 3 && type2 == 2 # YZ = iX
                phase *= 1.0im
            elseif type1 == 3 && type2 == 1 # YX = -iZ
                phase *= -1.0im
            elseif type1 == 1 && type2 == 3 # XY = iZ
                phase *= 1.0im
            end
        end
        active_bits >>= 1
        idx += 1
    end
    return (phase, PauliString(new_x, new_z))
end

# Conmutador [A, B]
function commutator(p1::PauliString, p2::PauliString)
    anticommutes = isodd(count_ones((p1.x_mask & p2.z_mask) ⊻ (p1.z_mask & p2.x_mask)))
    if anticommutes
        (phase_ab, p_ab) = multiply_paulis(p1, p2)
        return (2.0 * phase_ab, p_ab)
    else
        return (0.0im, PauliString(0,0))
    end
end