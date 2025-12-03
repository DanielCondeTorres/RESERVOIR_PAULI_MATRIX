"""
    add_shot_noise(value::Float64, n_shots::Int)

Simula el error estadístico de medir un observable con un número finito de tiros (shots).
Añade ruido gaussiano N(0, sigma) donde sigma depende de N_meas.
"""
function add_shot_noise(value::Float64, n_shots::Int)
    if n_shots <= 0; return value; end
    
    # Desviación estándar teórica (aproximada para Pauli Z)
    # sigma = sqrt((1 - value^2) / n_shots)
    # O la fórmula (11) si g es pequeño.
    sigma = sqrt(1.0 / n_shots) 
    
    noise = randn() * sigma
    return value + noise
end