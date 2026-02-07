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

function inject_state_EraseWrite(rho::Operator, qubit_idx::Int, rz::Float64; rx::Float64=0.0, ry::Float64=0.0)
    new_rho = Operator()
    bit_loc = 1 << qubit_idx
    tol = 1e-15
    has_x = abs(rx) > tol; has_y = abs(ry) > tol; has_z = abs(rz) > tol
    found_global_identity = false
    identity_mask = PauliString(0, 0)

    for (p, coeff) in rho
        if (p.x_mask >> qubit_idx) & 1 == 0 && (p.z_mask >> qubit_idx) & 1 == 0 # Check Pauli type 0 (Identity) manually
            if p == identity_mask; found_global_identity = true; end
            new_rho[p] = get(new_rho, p, 0.0im) + coeff
            if has_z; p_z = PauliString(p.x_mask, p.z_mask | bit_loc); new_rho[p_z] = get(new_rho, p_z, 0.0im) + (coeff * rz); end
            if has_x; p_x = PauliString(p.x_mask | bit_loc, p.z_mask); new_rho[p_x] = get(new_rho, p_x, 0.0im) + (coeff * rx); end
            if has_y; p_y = PauliString(p.x_mask | bit_loc, p.z_mask | bit_loc); new_rho[p_y] = get(new_rho, p_y, 0.0im) + (coeff * ry); end
        end
    end
    # Bloque de seguridad por si no había identidad
    if !found_global_identity
        if has_z; p_z = PauliString(0, bit_loc); new_rho[p_z] = get(new_rho, p_z, 0.0im) + rz; end
        if has_x; p_x = PauliString(bit_loc, 0); new_rho[p_x] = get(new_rho, p_x, 0.0im) + rx; end
        if has_y; p_y = PauliString(bit_loc, bit_loc); new_rho[p_y] = get(new_rho, p_y, 0.0im) + ry; end
    end
    return new_rho
end















function inject_state_EraseWrite2(rho::Operator, qubit_idx::Int, rz::Float64; rx::Float64=0.0, ry::Float64=0.0)
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









function inject_state_EraseWrite_matrix(rho::Matrix{ComplexF64}, qubit_idx::Int, 
    rz::Float64; rx::Float64=0.0, ry::Float64=0.0)
dim = size(rho, 1)
N_qubits = Int(log2(dim))

# 1. Matrices básicas
I2 = [1.0+0im 0.0; 0.0 1.0]
Zero_Proj = [1.0+0im 0.0; 0.0 0.0] # |0><0|
Flip_Op   = [0.0+0im 1.0; 0.0 0.0] # |0><1| (Pone 1 en 0)

# 2. Construir Operadores de Kraus para Resetear qubit k a |0>
# M0 = I...|0><0|...I
# M1 = I...|0><1|...I
# (Usamos qubit_idx + 1 porque Julia indexa desde 1)
k = qubit_idx + 1

M0 = (k==1) ? Zero_Proj : I2
M1 = (k==1) ? Flip_Op   : I2

for i in 2:N_qubits
M0 = kron(M0, (i==k) ? Zero_Proj : I2)
M1 = kron(M1, (i==k) ? Flip_Op   : I2)
end

# Aplicamos el Reset: rho -> |0>_k<0| (x) Tr_k(rho)
rho_reset = M0 * rho * M0' + M1 * rho * M1'

# 3. Rotar al estado objetivo (rx, ry, rz)
# Calculamos la unitaria que lleva |0> -> (rx, ry, rz)
norm = sqrt(rx^2 + ry^2 + rz^2)
if norm < 1e-9; return rho_reset; end # Si vector nulo, dejamos en |0>

# Ángulos esféricos
nx, ny, nz = rx/norm, ry/norm, rz/norm
theta = acos(nz)       # Polar
phi = atan(ny, nx)     # Azimutal

# U_rot = Rz(phi) * Ry(theta)
Rz = [exp(-im*phi/2) 0; 0 exp(im*phi/2)]
Ry = [cos(theta/2) -sin(theta/2); sin(theta/2) cos(theta/2)]
U_qubit = Rz * Ry

# Expandimos a todo el sistema
U_tot = (k==1) ? U_qubit : I2
for i in 2:N_qubits
U_tot = kron(U_tot, (i==k) ? U_qubit : I2)
end

return U_tot * rho_reset * U_tot'
end