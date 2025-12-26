using LinearAlgebra

# Asume que PauliString y Operator están definidos

"""
    build_parity_operator(n_qubits::Int) -> PauliString

Construye el operador de Paridad Total definido en las notas manuscritas (Pág 2).
Fórmula: P = Z_1 ⊗ Z_2 ⊗ ... ⊗ Z_N 

En representación de bits:
- x_mask = 00...0 (Ninguna X)
- z_mask = 11...1 (Z en todos los sitios)
"""
function build_parity_operator(n_qubits::Int)
    # Crear una máscara con N unos: (2^N) - 1
    # Ejemplo N=3: (1 << 3) - 1 = 1000 - 1 = 0111 (binario)
    z_mask_all = (1 << n_qubits) - 1
    
    return PauliString(0, z_mask_all)
end

"""
    measure_parity(rho::Operator, n_qubits::Int) -> Float64

Calcula el valor esperado de la paridad en el estado actual: ⟨P⟩ = Tr(ρ P).
- Si devuelve +1.0: El estado es totalmente par.
- Si devuelve -1.0: El estado es totalmente impar.
- Si devuelve 0.0: El estado es una mezcla o superposición sin paridad definida.
"""
function measure_parity(rho::Operator, n_qubits::Int)
    P = build_parity_operator(n_qubits)
    
    # El valor esperado es simplemente el coeficiente de P en la matriz de densidad
    # (o 0.0 si P no está explícitamente en el diccionario rho)
    parity_val = real(get(rho, P, 0.0im))
    
    return parity_val
end

"""
    check_symmetry_conservation(H::Operator, n_qubits::Int)

Verifica si el Hamiltoniano H conmuta con el operador de Paridad P.
Si [H, P] = 0, entonces la paridad se conserva y la dimensión efectiva se reduce a la mitad.
"""
function check_symmetry_conservation(H::Operator, n_qubits::Int)
    P = build_parity_operator(n_qubits)
    commutes = true
    
    # Verificamos término a término del Hamiltoniano
    for (h_term, coeff) in H
        # multiply_paulis devuelve la fase y el nuevo operador
        # Commutador: [A, B] = AB - BA
        # Dos Paulis conmutan si la fase de AB es igual a la de BA
        
        # Calculamos A * P
        (phase_ap, prod_ap) = multiply_paulis(h_term, P)
        
        # Calculamos P * A
        (phase_pa, prod_pa) = multiply_paulis(P, h_term)
        
        # Si AB != BA (fases opuestas), entonces anticonmutan
        if phase_ap != phase_pa
            println("⚠️ Violación de simetría encontrada en el término: $h_term")
            commutes = false
            # No hacemos break para ver todos los términos que violan
        end
    end
    
    if commutes
        println("✅ El Hamiltoniano CONSERVA la paridad P = Π Z_i.")
        println("   Nota: La dinámica ocurre en un subespacio de dimensión d/2.")
    else
        println("❌ El Hamiltoniano NO conserva la paridad P = Π Z_i.")
        println("   La paridad cambiará durante la evolución.")
    end
    
    return commutes
end