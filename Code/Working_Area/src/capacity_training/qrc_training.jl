using Statistics
using LinearAlgebra

"""
    train_reservoir(observables::Matrix{Float64}, targets::Vector{Float64}, washout::Int=20)

Entrena los pesos lineales W para ajustar: Target ≈ W * Observables.
Usa la ecuación (C2) del texto: ỹ_k = sum(w_m * O_k,m) + w_{L+1}
Implementa regresión lineal (Mínimos Cuadrados).

- observables: Matriz [T x L] (T pasos de tiempo, L observables medidos).
- targets: Vector [T] con los valores objetivo (ej: input retrasado).
- washout: Número de pasos iniciales a descartar (para olvidar condición inicial).

Retorna: (weights, predictions)
"""
function train_reservoir(observables::Matrix{Float64}, targets::Vector{Float64}, washout::Int=20)
    # 1. Descartar Washout (Transient)
    X = observables[washout+1:end, :]
    y = targets[washout+1:end]
    
    # 2. Añadir Bias (w_{L+1} en la ecuación)
    # Añadimos una columna de 1s a la matriz de observables
    rows = size(X, 1)
    X_bias = hcat(X, ones(rows))
    
    # 3. Resolver Mínimos Cuadrados: w = (X'X)^-1 X'y
    # El operador \ de Julia hace esto eficientemente y estable.
    weights = X_bias \ y
    
    # 4. Generar Predicciones (In-Sample o Test según uso)
    predictions = X_bias * weights
    
    return weights, predictions
end

"""
    calculate_capacity(target::Vector{Float64}, prediction::Vector{Float64})

Calcula la Capacidad de Memoria (C) según la Ec. (C3):
C = cov(y, ỹ)^2 / (var(y) * var(ỹ))
Esto es equivalente al Coeficiente de Determinación R^2 o correlación de Pearson al cuadrado.
"""
function calculate_capacity(target::Vector{Float64}, prediction::Vector{Float64})
    # Asegurar que tengan la misma longitud
    n = min(length(target), length(prediction))
    y = target[end-n+1:end]
    y_pred = prediction[end-n+1:end]
    
    # Calcular covarianza y varianzas
    cv = cov(y, y_pred)
    v_y = var(y)
    v_pred = var(y_pred)
    
    # Evitar división por cero si es constante
    if v_y < 1e-12 || v_pred < 1e-12
        return 0.0
    end
    
    # Fórmula (C3)
    C = (cv^2) / (v_y * v_pred)
    return C
end