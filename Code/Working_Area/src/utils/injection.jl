# ==========================================
# INYECCIÓN DE INPUT (SCHRÖDINGER / FORWARD)
# ==========================================

function get_pauli_type(p::PauliString, qubit_idx::Int)
    x_bit = (p.x_mask >> qubit_idx) & 1
    z_bit = (p.z_mask >> qubit_idx) & 1
    return x_bit + 2 * z_bit 
end

"""
    inject_input_schrodinger(rho::Operator, qubit_idx::Int, u::Float64)
    
Versión Corregida: Maneja correctamente el caso donde no hay Identidad explícita.
"""

function apply_global_dephasing(rho::Operator, g::Float64)
    new_rho = Operator()
    
    # Factor de decaimiento: (1 - 2g) o similar dependiendo de la definición exacta del canal.
    # Asumimos canal de despolarización de fase estándar: rho -> (1-g)rho + g X rho X
    # Esto deja X quieto, y multiplica Y y Z por (1-2g).
    decay_factor = (1.0 - 2.0 * g)
    
    for (p, coeff) in rho
        # Calculamos cuántos Y o Z hay en total en esta cadena
        # (Son los bits encendidos en z_mask, ya que Z=01 y Y=11, ambos tienen el bit z activado)
        n_non_commuting = count_ones(p.z_mask)
        
        # Cada operador no-conmutativo aporta un factor de decaimiento
        total_decay = decay_factor ^ n_non_commuting
        
        if abs(coeff * total_decay) > 1e-20
            new_rho[p] = get(new_rho, p, 0.0im) + (coeff * total_decay)
        end
    end
    
    return new_rho
end



# ==============================================================================
#  FUNCIONES DE INYECCIÓN DE ESTADOS (Reset & Write)
# ==============================================================================

"""
    get_pauli_type(p::PauliString, qubit_idx::Int) -> Int

Helper para identificar qué operador hay en un qubit específico:
0 = I, 1 = X, 2 = Z, 3 = Y
"""
function get_pauli_type(p::PauliString, qubit_idx::Int)
    x_bit = (p.x_mask >> qubit_idx) & 1
    z_bit = (p.z_mask >> qubit_idx) & 1
    return x_bit + 2 * z_bit 
end

"""
    inject_state(rho::Operator, qubit_idx::Int, rz::Float64; rx::Float64=0.0, ry::Float64=0.0)

Inyecta un estado en el qubit `qubit_idx` definido por el vector de Bloch (rx, ry, rz).
Por defecto, asume que el input principal `rz` es el dato `u` y que `rx` y `ry` son 0.

Uso típico (Reservoir Computing):
    inject_state(rho, 1, u)  -> Inyecta u en Z, nada en X/Y.

Uso avanzado (Ruido o Control):
    inject_state(rho, 1, u, rx=0.1) -> Inyecta u en Z y 0.1 en X.
"""
function inject_state(rho::Operator, qubit_idx::Int, rz::Float64; rx::Float64=0.0, ry::Float64=0.0)
    new_rho = Operator()
    bit_loc = 1 << qubit_idx
    tol = 1e-15
    
    # Pre-calculamos si necesitamos procesar X o Y para ahorrar tiempo (Optimización)
    has_x = abs(rx) > tol
    has_y = abs(ry) > tol
    has_z = abs(rz) > tol

    # 1. Procesar términos existentes (Correlaciones)
    for (p, coeff) in rho
        # Chequeo rápido: si el qubit no es identidad, el término muere (Traza parcial)
        # get_pauli_type: 0=I, 1=X, 2=Z, 3=Y
        if get_pauli_type(p, qubit_idx) == 0
            
            # A) Mantenemos la memoria (término con I en este qubit)
            new_rho[p] = get(new_rho, p, 0.0im) + coeff
            
            # B) Inyectamos Z (El input estándar)
            if has_z
                p_z = PauliString(p.x_mask, p.z_mask | bit_loc)
                new_rho[p_z] = get(new_rho, p_z, 0.0im) + (coeff * rz)
            end

            # C) Inyectamos X (Solo si rx != 0)
            if has_x
                p_x = PauliString(p.x_mask | bit_loc, p.z_mask)
                new_rho[p_x] = get(new_rho, p_x, 0.0im) + (coeff * rx)
            end

            # D) Inyectamos Y (Solo si ry != 0)
            if has_y
                p_y = PauliString(p.x_mask | bit_loc, p.z_mask | bit_loc)
                new_rho[p_y] = get(new_rho, p_y, 0.0im) + (coeff * ry)
            end
        end
    end
    
    # 2. Inyectar el término fresco (donde el resto de la cadena era Identidad implícita)
    if has_z
        p_fresh_z = PauliString(0, bit_loc)
        new_rho[p_fresh_z] = get(new_rho, p_fresh_z, 0.0im) + rz
    end
    
    if has_x
        p_fresh_x = PauliString(bit_loc, 0)
        new_rho[p_fresh_x] = get(new_rho, p_fresh_x, 0.0im) + rx
    end
    
    if has_y
        p_fresh_y = PauliString(bit_loc, bit_loc)
        new_rho[p_fresh_y] = get(new_rho, p_fresh_y, 0.0im) + ry
    end
    
    return new_rho
end