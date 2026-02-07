# src/utils/measurements.jl
using LinearAlgebra

"""
    build_paper_basis(n_qubits::Int) -> Vector{PauliString}
Construye la lista exacta de observables del paper: Z_i, X_i, y Z_i Z_{i+1}.
"""
function build_paper_basis(n_qubits::Int)
    # 1. Z locales + 2. X locales + 3. Correlaciones ZZ vecinos
    return vcat(
        [PauliString(0, 1 << i) for i in 0:(n_qubits-1)],          # Z_i
        [PauliString(1 << i, 0) for i in 0:(n_qubits-1)],          # X_i
        [PauliString(0, (1<<i)|(1<<(i+1))) for i in 0:(n_qubits-2)] # Z_i Z_{i+1}
    )
end
"""
    measure_observable(rho, P; projective=false) -> (val, rho_new)
Mide un observable P.
- projective=false: Retorna ⟨P⟩ y el rho intacto.
- projective=true:  Retorna ±1 y el rho colapsado (Backaction).
"""
function measure_observable(rho::Operator, P::PauliString; projective::Bool=false)
    # 1. Valor esperado ⟨P⟩ = Tr(P ρ)
    expectation = real(get(rho, P, 0.0im))
    if !projective
        return expectation, rho
    end
    # 2. Lógica Proyectiva (Colapso)
    prob_plus = clamp((1.0 + expectation) / 2.0, 0.0, 1.0)
    outcome = rand() < prob_plus ? 1.0 : -1.0
    
    # Backaction: ρ' = Π ρ Π / p
    # Simplificado: ρ' = (ρ + o·{P, ρ}) / (2p) donde {P,ρ} = Pρ + ρP
    # Nota: Implementación completa requeriría lógica de multiplicación de Paulis similar a la tuya anterior.
    # Si solo vas a usar Expectation Values para la Fig 3a, esto es suficiente.
    # Para no complicar, aquí asumimos que NO necesitas backaction completa para la Fig 3a.
    return outcome, rho # (Devolvemos rho sin colapsar si no implementamos apply_backaction completo aquí)
end

"""
    extract_features(rho, observables) -> Vector{Float64}
Mide una lista de observables y devuelve el vector de resultados.
"""
function extract_features(rho::Operator, observables::Vector{PauliString})
    # Mapea la función measure_observable sobre la lista y extrae solo el valor (índice 1)
    return [measure_observable(rho, P; projective=false)[1] for P in observables]
end













"""
    measure_observable_matrix(rho, op; projective::Bool=false, gamma::Float64=1.0)

Realiza una medición del operador 'op' sobre el estado 'rho'.

Argumentos:
- `rho`: Matriz de densidad actual.
- `op`: Observable a medir (matriz hermítica con autovalores ±1).
- `projective`: 
    - false: Devuelve valor esperado <op> y NO altera rho (Ideal para Machine Learning).
    - true:  Simula un colapso físico (Backaction).
- `gamma`: Fuerza de la medida (solo si projective=true). 
    - 1.0: Medida Proyectiva Fuerte (Colapso total a |+> o |->).
    - 0.0: Medida Nula (No hace nada, resultado aleatorio).
    - 0.5: Medida Débil (Colapso parcial).

Retorna:
- (outcome, rho_new)
"""
function measure_observable_matrix(rho::Matrix{ComplexF64}, 
    op::Matrix{ComplexF64}, 
    gamma::Float64=1.0; # Ahora es posicional con default 1.0
    projective::Bool=false) # Keyword argument    
    # 1. MODO LECTURA (Machine Learning / ESP)
    # Calculamos el valor esperado <op> = Tr(rho * op)
    expectation = real(tr(rho * op))
    
    if !projective
        # Devolvemos el valor continuo y el estado intacto
        return expectation, rho 
    end

    # ---------------------------------------------------------
    # 2. MODO COLAPSO FÍSICO (Generalizado con Kraus)
    # ---------------------------------------------------------
    
    # Dimensiones e Identidad
    dim = size(rho, 1)
    I_mat = Matrix{ComplexF64}(I, dim, dim)
    
    # Definimos los Proyectores Ortogonales (Tus notas Π±)
    P_plus  = (I_mat + op) / 2.0
    P_minus = (I_mat - op) / 2.0
    
    # Coeficientes para Operadores de Kraus Débiles (Generalización)
    # M± = sqrt((1±γ)/2) Π+ + sqrt((1∓γ)/2) Π-
    # Si gamma=1: g_main=1, g_cross=0 -> M+ = P+, M- = P- (Proyectiva pura)
    g_main  = sqrt((1.0 + gamma) / 2.0)
    g_cross = sqrt((1.0 - gamma) / 2.0)
    
    # Construimos los Operadores de Kraus para los dos posibles resultados
    M_plus  = (g_main * P_plus) + (g_cross * P_minus)
    M_minus = (g_cross * P_plus) + (g_main * P_minus)
    
    # Calculamos la probabilidad REAL de obtener el resultado +1
    # P(+1) = Tr(M+ rho M+†)
    # Nota: Si gamma < 1, la probabilidad se acerca a 0.5 (más ruido)
    rho_unnorm_plus = M_plus * rho * M_plus'
    prob_plus_outcome = real(tr(rho_unnorm_plus))
    
    # --- Tirar la moneda (Monte Carlo) ---
    if rand() < prob_plus_outcome
        outcome = 1.0
        # El estado colapsa hacia +1 (con fuerza gamma)
        # Normalizamos: rho' = (M+ rho M+') / P(+)
        rho_new = rho_unnorm_plus / prob_plus_outcome
    else
        outcome = -1.0
        # El estado colapsa hacia -1 (con fuerza gamma)
        rho_unnorm_minus = M_minus * rho * M_minus'
        prob_minus_outcome = real(tr(rho_unnorm_minus))
        
        # Evitar división por cero en casos extremos
        if prob_minus_outcome < 1e-12
            rho_new = rho # Fallback raro
        else
            rho_new = rho_unnorm_minus / prob_minus_outcome
        end
    end
    
    return outcome, rho_new
end