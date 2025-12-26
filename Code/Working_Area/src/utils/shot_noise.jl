# src/utils/shot_noise.jl
using Random

"""
    apply_shot_noise(val::Float64, n_shots::Real) -> Float64

Aplica ruido de disparo a un valor esperado individual.
S=1.5e6 es el valor usado por Nathan para observables locales.
"""
function apply_shot_noise(val::Float64, n_shots::Real)
    if n_shots <= 0
        return val
    end
    # La varianza de la media de S disparos es (1 - <P>^2) / S
    # Usamos max(0.0, ...) para evitar errores numéricos de precisión
    sigma = sqrt(max(0.0, 1.0 - val^2) / n_shots)
    # Añadimos ruido Gaussiano
    noisy_val = val + randn() * sigma
    # Aseguramos que el valor físico permanezca en el rango [-1, 1]
    return clamp(noisy_val, -1.0, 1.0)
end

"""
    add_shot_noise_to_series(series::Vector{Float64}, n_shots::Real)
    
Devuelve una nueva serie con ruido aplicado a cada punto temporal.
"""
function add_shot_noise_to_series(series::Vector{Float64}, n_shots::Real)
    return [apply_shot_noise(v, n_shots) for v in series]
end