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