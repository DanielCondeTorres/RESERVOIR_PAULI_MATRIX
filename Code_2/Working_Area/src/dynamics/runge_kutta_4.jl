# --- FUNCIÓN DE APOYO: Derivada (Von Neumann) ---
function derivative(O::Operator, H::Operator)
    # dO/dt = -i [H, O]
    res = Operator()
    for (ph, ch) in H, (po, co) in O
        p_ab, f_ab = multiply_paulis(ph, po)
        p_ba, f_ba = multiply_paulis(po, ph)
        val = -im * ((ch * co * f_ab) - (ch * co * f_ba))
        if abs(val) > 1e-15
            res[p_ab] = get(res, p_ab, 0.0im) + val
        end
    end
    return res
end

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