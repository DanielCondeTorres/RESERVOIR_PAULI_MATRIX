# ==============================================================================
# src/utils/pauli_algebra.jl
# Operaciones fundamentales de álgebra lineal sobre PauliStrings y Operadores
# ==============================================================================

const TOLERANCE = 1e-15

"""
    multiply_paulis(p1, p2) -> (phase, p_new)
Producto de dos cadenas de Pauli.
"""
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
            
            if t1==1 && t2==2; phase *= -1.0im
            elseif t1==2 && t2==1; phase *= 1.0im
            elseif t1==2 && t2==3; phase *= -1.0im
            elseif t1==3 && t2==2; phase *= 1.0im
            elseif t1==3 && t2==1; phase *= -1.0im
            elseif t1==1 && t2==3; phase *= 1.0im
            end
        end
        active_bits >>= 1
        idx += 1
    end
    return (phase, PauliString(new_x, new_z))
end

"""
    commutator(p1, p2) -> (val, p_new)
Calcula [P1, P2]. Retorna 0 si conmutan.
"""
function commutator(p1::PauliString, p2::PauliString)
    if isodd(count_ones((p1.x_mask & p2.z_mask) ⊻ (p1.z_mask & p2.x_mask)))
        (ph, pn) = multiply_paulis(p1, p2)
        return (2.0 * ph, pn)
    end
    return (0.0im, PauliString(0,0))
end

"""
    add_ops(A, B, scale_B)
Suma de operadores: C = A + (scale * B)
"""
function add_ops(A::Operator, B::Operator, scale_B::ComplexF64)
    C = copy(A)
    for (p, val) in B
        new_val = get(C, p, 0.0im) + val * scale_B
        if abs(new_val) > TOLERANCE 
            C[p] = new_val
        else
            delete!(C, p)
        end
    end
    return C
end


#function truncate_operator!(O::Operator, max_terms::Int)
#    if length(O) <= max_terms; return; end
#    all_terms = collect(O)
#    sort!(all_terms, by = x -> abs(x[2]), rev = true)
#for i in 1:max_terms
    #    empty!(O)
    #    O[all_terms[i][1]] = all_terms[i][2]
    #end
#end





"""
    truncate_operator!(O::Operator, max_terms::Int)

Implementa la estrategia de truncamiento descrita en la sección 'Simulation complexity' del texto:
"we will truncate after the largest M Pauli coefficients."

Utiliza `partialsort!` para encontrar los Top-M términos sin ordenar todo el vector (O(N) vs O(N log N)),
lo cual es crucial para la eficiencia cuando el espacio de Hilbert es grande (N=20).
"""
function truncate_operator!(O::Operator, max_terms::Int=typemax(Int))
    # 1. Si el operador ya es pequeño, no hacemos nada (ahorro de tiempo)
    if length(O) <= max_terms
        return
    end
    
    # 2. Convertir el diccionario a un vector de pares (PauliString, Coeff)
    all_terms = collect(O)
    
    # 3. SELECCIÓN DE LOS TOP-M (Optimizada según el texto)
    # Usamos partialsort! para encontrar solo los 'max_terms' coeficientes más grandes
    # en valor absoluto (abs(x[2])), ordenados de mayor a menor (rev=true).
    # Esto evita ordenar la "cola" de coeficientes pequeños que vamos a borrar.
    top_terms = partialsort!(all_terms, 1:max_terms, by = x -> abs(x[2]), rev = true)
    
    # 4. RECONSTRUCCIÓN
    # Vaciamos el operador original y reinsertamos solo los supervivientes
    empty!(O)
    for (p, c) in top_terms
        O[p] = c
    end
end