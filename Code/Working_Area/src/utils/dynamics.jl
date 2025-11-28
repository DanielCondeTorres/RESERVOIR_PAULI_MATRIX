# ==========================================
# UTILIDADES DINÁMICAS (INTEGRADOR)
# ==========================================

# Suma de operadores con truncamiento (limpieza de memoria)
function add_ops(A::Operator, B::Operator, scale_B::ComplexF64)
    C = copy(A)
    for (p, val) in B
        new_val = get(C, p, 0.0im) + val * scale_B
        if abs(new_val) > 1e-10 
            C[p] = new_val
        else
            delete!(C, p)
        end
    end
    return C
end

# Ecuación de Heisenberg: dO/dt = i [H, O]
function derivative(O::Operator, H::Operator)
    dO = Operator()
    for (p_op, c_op) in O
        for (p_H, c_H) in H
            (val_comm, p_new) = commutator(p_H, p_op)
            if val_comm != 0
                coeff = im * c_H * c_op * val_comm
                if abs(coeff) > 1e-15
                    dO[p_new] = get(dO, p_new, 0.0im) + coeff
                end
            end
        end
    end
    return dO
end

# Solver Runge-Kutta 4 (RK4)
function step_rk4(O::Operator, H::Operator, dt::Float64)
    k1 = derivative(O, H)
    k2 = derivative(add_ops(O, k1, 0.5*dt + 0.0im), H)
    k3 = derivative(add_ops(O, k2, 0.5*dt + 0.0im), H)
    k4 = derivative(add_ops(O, k3, dt + 0.0im), H)
    
    O_new = add_ops(O, k1, dt/6.0 + 0.0im)
    O_new = add_ops(O_new, k2, dt/3.0 + 0.0im)
    O_new = add_ops(O_new, k3, dt/3.0 + 0.0im)
    O_new = add_ops(O_new, k4, dt/6.0 + 0.0im)
    return O_new
end



"""
    truncate_operator!(O::Operator, max_terms::Int)

Reduce el operador 'O' para quedarse solo con los 'max_terms' términos 
que tienen los coeficientes más grandes (en valor absoluto).
Esto simula la aproximación de "Pauli Truncation".
"""
function truncate_operator!(O::Operator, max_terms::Int)
    # Si tenemos menos términos que el límite, no hacemos nada
    if length(O) <= max_terms
        return O
    end
    
    # 1. Extraer todos los pares (Pauli => Coeficiente)
    all_terms = collect(O)
    
    # 2. Ordenar por magnitud del coeficiente (de mayor a menor)
    #    Usamos abs(x[2]) porque x[2] es el coeficiente complejo
    sort!(all_terms, by = x -> abs(x[2]), rev = true)
    
    # 3. Limpiar el operador original y rellenar solo con los Top-K
    empty!(O)
    for i in 1:max_terms
        (p, c) = all_terms[i]
        O[p] = c
    end
    
    return O
end













# src/utils/dynamics.jl

# Multiplicación Pauli
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
            # Tabla simplificada
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

function commutator(p1::PauliString, p2::PauliString)
    if isodd(count_ones((p1.x_mask & p2.z_mask) ⊻ (p1.z_mask & p2.x_mask)))
        (ph, pn) = multiply_paulis(p1, p2)
        return (2.0 * ph, pn)
    end
    return (0.0im, PauliString(0,0))
end

function add_ops(A::Operator, B::Operator, scale::ComplexF64)
    C = copy(A)
    for (p, val) in B
        new_val = get(C, p, 0.0im) + val * scale
        if abs(new_val) > 1e-20; C[p] = new_val; else; delete!(C, p); end
    end
    return C
end

function derivative(O::Operator, H::Operator)
    dO = Operator()
    for (pop, cop) in O, (ph, ch) in H
        (val, pnew) = commutator(ph, pop)
        if val != 0
            c_new = im * ch * cop * val
            if abs(c_new) > 1e-20; dO[pnew] = get(dO, pnew, 0.0im) + c_new; end
        end
    end
    return dO
end

function step_rk4(O, H, dt)
    k1 = derivative(O, H)
    k2 = derivative(add_ops(O, k1, 0.5*dt + 0.0im), H)
    k3 = derivative(add_ops(O, k2, 0.5*dt + 0.0im), H)
    k4 = derivative(add_ops(O, k3, dt + 0.0im), H)
    
    O = add_ops(O, k1, dt/6.0 + 0.0im)
    O = add_ops(O, k2, dt/3.0 + 0.0im)
    O = add_ops(O, k3, dt/3.0 + 0.0im)
    O = add_ops(O, k4, dt/6.0 + 0.0im)
    return O
end

"""
    truncate_operator!(O, max_terms)
    Recorta el operador manteniendo solo los términos más grandes.
"""
function truncate_operator!(O::Operator, max_terms::Int)
    if length(O) <= max_terms; return; end
    all = collect(O)
    sort!(all, by=x->abs(x[2]), rev=true)
    empty!(O)
    for i in 1:max_terms
        O[all[i][1]] = all[i][2]
    end
end