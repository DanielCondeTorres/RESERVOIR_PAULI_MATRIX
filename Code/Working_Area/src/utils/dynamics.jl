# ==============================================================================
# src/utils/dynamics.jl
# Integradores numéricos y ecuación de Heisenberg
# Requiere: pauli_algebra.jl cargado previamente
# ==============================================================================

"""
    derivative(O, H)
Calcula dO/dt = i [H, O] (Ecuación de Heisenberg)
"""
function derivative(O::Operator, H::Operator)
    dO = Operator()
    for (pop, cop) in O, (ph, ch) in H
        # Llama a commutator definido en pauli_algebra.jl
        (val_comm, pnew) = commutator(ph, pop)
        
        if val_comm != 0
            c_new = im * ch * cop * val_comm
            if abs(c_new) > TOLERANCE
                dO[pnew] = get(dO, pnew, 0.0im) + c_new
            end
        end
    end
    return dO
end

"""
    step_rk4(O, H, dt)
Integración temporal Runge-Kutta 4.
"""
function step_rk4(O::Operator, H::Operator, dt::Float64)
    # Llama a add_ops y derivative
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