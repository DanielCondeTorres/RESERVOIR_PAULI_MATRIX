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
function inject_input_schrodinger(rho::Operator, qubit_idx::Int, u::Float64)
    new_rho = Operator()
    
    # 1. Procesar términos existentes (Correlaciones pasadas)
    for (p, coeff) in rho
        type = get_pauli_type(p, qubit_idx)
        
        # Solo sobreviven los términos que NO tienen X, Y o Z en el qubit que estamos tocando.
        # (Si tiene Z, al trazarlo da 0, así que se elimina. Si tiene I, sobrevive).
        if type == 0 # Identidad (I) en este qubit
            
            # A) Mantener la correlación antigua (El "recuerdo" en los otros qubits)
            new_rho[p] = get(new_rho, p, 0.0im) + coeff
            
            # B) Mezclar el nuevo input con ese recuerdo (No-linealidad)
            # Añadimos u * Z_qubit * (Resto_de_la_cadena)
            z_mask_bit = (1 << qubit_idx)
            p_con_z = PauliString(p.x_mask, p.z_mask | z_mask_bit)
            
            new_rho[p_con_z] = get(new_rho, p_con_z, 0.0im) + (coeff * u)
        end
    end
    
    # 2. INYECCIÓN PURA (El término nuevo fresco)
    # Como el diccionario no guarda la "Identidad Global" explícitamente,
    # tenemos que añadir manualmente el input u * Z_0.
    
    p_fresh = PauliString(0, 1 << qubit_idx) # Z en el qubit indicado
    new_rho[p_fresh] = get(new_rho, p_fresh, 0.0im) + u
    
    return new_rho
end




# ... (Mantén la función get_pauli_type anterior) ...

"""
    inject_general_bloch_state(rho::Operator, qubit_idx::Int, rx::Float64, ry::Float64, rz::Float64)

Implementa el mapa Erase & Write genérico:
1. Borra el qubit_idx (Traza parcial).
2. Escribe el estado mixto definido por el vector de Bloch (rx, ry, rz).
   rho_new = rho_rest ⊗ (I + rx*X + ry*Y + rz*Z) / 2
"""
function inject_general_bloch_state(rho::Operator, qubit_idx::Int, rx::Float64, ry::Float64, rz::Float64)
    new_rho = Operator()
    
    # Máscaras de bits para operadores X, Z (y Y es X|Z)
    bit_loc = 1 << qubit_idx
    
    for (p, coeff) in rho
        # Verificamos qué hay en el qubit que vamos a tocar
        type = get_pauli_type(p, qubit_idx)
        
        # Solo sobreviven los términos que tienen Identidad (I) en este qubit.
        # (Esto equivale a borrar la información vieja de este qubit).
        if type == 0 
            # 1. Parte I: Se mantiene el resto de la cadena (rho_rest ⊗ I)
            new_rho[p] = get(new_rho, p, 0.0im) + coeff
            
            # 2. Parte X: rho_rest ⊗ X
            if abs(rx) > 1e-15
                p_x = PauliString(p.x_mask | bit_loc, p.z_mask)
                new_rho[p_x] = get(new_rho, p_x, 0.0im) + (coeff * rx)
            end

            # 3. Parte Y: rho_rest ⊗ Y
            if abs(ry) > 1e-15
                p_y = PauliString(p.x_mask | bit_loc, p.z_mask | bit_loc)
                new_rho[p_y] = get(new_rho, p_y, 0.0im) + (coeff * ry)
            end

            # 4. Parte Z: rho_rest ⊗ Z
            if abs(rz) > 1e-15
                p_z = PauliString(p.x_mask, p.z_mask | bit_loc)
                new_rho[p_z] = get(new_rho, p_z, 0.0im) + (coeff * rz)
            end
        end
    end
    
    # 5. Inyectar el término fresco (por si rho estaba vacío)
    if abs(rx) > 1e-15
        p_fresh_x = PauliString(bit_loc, 0)
        new_rho[p_fresh_x] = get(new_rho, p_fresh_x, 0.0im) + rx
    end
    if abs(ry) > 1e-15
        p_fresh_y = PauliString(bit_loc, bit_loc)
        new_rho[p_fresh_y] = get(new_rho, p_fresh_y, 0.0im) + ry
    end
    if abs(rz) > 1e-15
        p_fresh_z = PauliString(0, bit_loc)
        new_rho[p_fresh_z] = get(new_rho, p_fresh_z, 0.0im) + rz
    end
    
    return new_rho
end



"""
    apply_global_dephasing(rho, strength_g)
    
    Aplica ruido de desfase en la base X a TODOS los qubits.
    Físicamente: Los términos que tienen Y o Z decaen por un factor.
    Esto simula la pérdida de coherencia o medición débil descrita en el paper.
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






# src/utils/injection.jl

function get_pauli_type(p::PauliString, qubit_idx::Int)
    x_bit = (p.x_mask >> qubit_idx) & 1
    z_bit = (p.z_mask >> qubit_idx) & 1
    return x_bit + 2 * z_bit 
end

function inject_general_bloch_state(rho::Operator, qubit_idx::Int, rx::Float64, ry::Float64, rz::Float64)
    new_rho = Operator()
    bit_loc = 1 << qubit_idx
    
    for (p, coeff) in rho
        type = get_pauli_type(p, qubit_idx)
        if type == 0 
            new_rho[p] = get(new_rho, p, 0.0im) + coeff
            if abs(rx) > 1e-15
                px = PauliString(p.x_mask | bit_loc, p.z_mask)
                new_rho[px] = get(new_rho, px, 0.0im) + (coeff * rx)
            end
            if abs(rz) > 1e-15
                pz = PauliString(p.x_mask, p.z_mask | bit_loc)
                new_rho[pz] = get(new_rho, pz, 0.0im) + (coeff * rz)
            end
        end
    end
    if abs(rx) > 1e-15; px = PauliString(bit_loc, 0); new_rho[px] = get(new_rho, px, 0.0im) + rx; end
    if abs(rz) > 1e-15; pz = PauliString(0, bit_loc); new_rho[pz] = get(new_rho, pz, 0.0im) + rz; end
    
    return new_rho
end

"""
    apply_global_dephasing(rho, g)
    Aplica ruido en base X. Los términos que anticonmutan con X (Y, Z) decaen.
"""
function apply_global_dephasing(rho::Operator, g::Float64)
    new_rho = Operator()
    decay_factor = (1.0 - 2.0 * g) # Factor de reducción para Y y Z
    
    for (p, coeff) in rho
        # Contamos cuántos Y o Z hay (bits z encendidos)
        n_anticommute = count_ones(p.z_mask)
        total_decay = decay_factor ^ n_anticommute
        
        new_val = coeff * total_decay
        if abs(new_val) > 1e-20
            new_rho[p] = new_val
        end
    end
    return new_rho
end