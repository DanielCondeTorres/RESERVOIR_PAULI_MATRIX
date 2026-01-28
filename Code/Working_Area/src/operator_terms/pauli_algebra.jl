# ==============================================================================
# src/operator_terms/pauli_algebra.jl
# Operaciones fundamentales de álgebra lineal sobre PauliStrings y Operadores
# ==============================================================================

const TOLERANCE = 1e-15

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

# ==============================================================================
# 1. TRUNCAMIENTO (¡RECUPERADA!)
# ==============================================================================
"""
    truncate_operator!(O::Operator, max_terms::Int)
    
Mantiene solo los 'max_terms' coeficientes más grandes del operador.
Esencial para evitar que la memoria explote en simulaciones largas.
"""
function truncate_operator!(O::Operator, max_terms::Int)
    # Si el operador es pequeño, no hacemos nada
    if length(O) <= max_terms
        return
    end
    
    # 1. Convertir a vector y ordenar por magnitud (mayor a menor)
    all_terms = collect(O)
    sort!(all_terms, by = x -> abs(x[2]), rev = true)
    
    # 2. Vaciar el diccionario
    empty!(O)
    
    # 3. Reinsertar solo los top-M términos
    for i in 1:max_terms
        (p, val) = all_terms[i]
        O[p] = val
    end
end

# ==============================================================================
# 2. ÁLGEBRA (Add, Multiply, Commutator)
# ==============================================================================

function add_ops!(target::Operator, source::Operator, scale::Number)
    if abs(scale) < TOLERANCE; return target; end

    for (p, val_source) in source
        val_current = get(target, p, 0.0im)
        new_val = val_current + (val_source * scale)
        
        if abs(new_val) > TOLERANCE 
            target[p] = new_val
        else
            delete!(target, p)
        end
    end
    return target
end

# Sobrecarga para crear copia
function add_ops(A::Operator, B::Operator, scale::Number)
    C = copy(A)
    add_ops!(C, B, scale)
    return C
end

function multiply_paulis(p1::PauliString, p2::PauliString)
    new_x = p1.x_mask ⊻ p2.x_mask
    new_z = p1.z_mask ⊻ p2.z_mask
    phase = 1.0 + 0.0im
    
    active_bits = (p1.x_mask | p1.z_mask) & (p2.x_mask | p2.z_mask)
    idx = 0
    while active_bits > 0
        if (active_bits & 1) == 1
            t1 = ((p1.x_mask >> idx) & 1) + 2 * ((p1.z_mask >> idx) & 1)
            t2 = ((p2.x_mask >> idx) & 1) + 2 * ((p2.z_mask >> idx) & 1)
            
            if t1 != t2
                if (t1==1 && t2==2) || (t1==2 && t2==3) || (t1==3 && t2==1) 
                    phase *= -1.0im
                else
                    phase *= 1.0im
                end
            end
        end
        active_bits >>= 1
        idx += 1
    end
    return (phase, PauliString(new_x, new_z))
end

function commutator(p1::PauliString, p2::PauliString)
    anticommutes = isodd(count_ones((p1.x_mask & p2.z_mask) ⊻ (p1.z_mask & p2.x_mask)))
    if anticommutes
        (ph, pn) = multiply_paulis(p1, p2)
        return (2.0 * ph, pn)
    end
    return (0.0im, PauliString(0,0))
end

# ==============================================================================
# 3. CREACIÓN DE OBSERVABLES
# ==============================================================================

function create_pauli_observable(tipo_str::String, lista_qubits::Vector{Int})
    total_x = 0
    total_z = 0
    
    tipo = uppercase(tipo_str)
    usar_x = (tipo == "X" || tipo == "Y")
    usar_z = (tipo == "Z" || tipo == "Y")
    
    for q in lista_qubits
        shift = q - 1 
        if usar_x; total_x |= (1 << shift); end
        if usar_z; total_z |= (1 << shift); end
    end
    
    return PauliString(total_x, total_z)
end





function operator_to_dense_matrix(H_dict::Operator, n_qubits::Int)
    dim = 2^n_qubits
    H_mat = zeros(ComplexF64, dim, dim)
    I2, X, Z, Y = [1 0; 0 1], [0 1; 1 0], [1 0; 0 -1], [0 -im; im 0]

    for (pauli, coeff) in H_dict
        # Empezamos con la matriz de identidad de 1x1
        term_mat = ComplexF64[1.0;;] 
        for k in 0:(n_qubits - 1)
            x = (pauli.x_mask >> k) & 1
            z = (pauli.z_mask >> k) & 1
            
            gate = I2
            if x == 1 && z == 0; gate = X
            elseif x == 0 && z == 1; gate = Z
            elseif x == 1 && z == 1; gate = Y
            end
            
            # kron(Qubits_posteriores, Qubits_anteriores)
            term_mat = kron(gate, term_mat)
        end
        H_mat += coeff * term_mat
    end
    return H_mat
end