"""
    separability_metrics.jl

Módulo para calcular métricas de separabilidad en Quantum Reservoir Computing.
Basado en:
- Fujii & Nakajima (2021): https://arxiv.org/abs/2102.11831
- Domingo et al. (2024): https://arxiv.org/pdf/2403.02686
"""

using LinearAlgebra

# ==============================================================================
# MÉTRICAS DE DISTANCIA ENTRE ESTADOS
# ==============================================================================

"""
    euclidean_distance(vec1::Vector, vec2::Vector)

Distancia euclidiana entre dos vectores de observables.
"""
function euclidean_distance(vec1::Vector, vec2::Vector)
    return sqrt(sum((vec1 .- vec2).^2))
end

"""
    trace_distance(rho1::Dict, rho2::Dict)

Distancia de traza entre dos estados cuánticos: D(ρ₁,ρ₂) = (1/2)||ρ₁ - ρ₂||₁
Aproximación usando los operadores disponibles.
"""
function trace_distance(rho1::Dict, rho2::Dict)
    # Para estados puros esto se simplifica considerablemente
    # Aquí calculamos la norma Frobenius como aproximación
    dist = 0.0
    all_keys = union(keys(rho1), keys(rho2))
    
    for key in all_keys
        val1 = get(rho1, key, 0.0)
        val2 = get(rho2, key, 0.0)
        dist += abs(val1 - val2)^2
    end
    
    return sqrt(dist)
end

# ==============================================================================
# SEPARABILIDAD INSTANTÁNEA
# ==============================================================================

"""
    instantaneous_separation(obs_current::Vector, obs_previous::Vector)

Separación instantánea entre estados consecutivos.
Se calcula cuando hay transición de input.
"""
function instantaneous_separation(obs_current::Vector, obs_previous::Vector)
    return euclidean_distance(obs_current, obs_previous)
end

# ==============================================================================
# KERNEL DE SEPARABILIDAD (SP)
# ==============================================================================

"""
    separation_kernel(trajectory::Matrix, inputs::Vector; window::Int=1)

Calcula el kernel de separabilidad según Fujii & Nakajima.
SP(u,v) mide cuán separables son dos secuencias de inputs u y v.

# Argumentos
- `trajectory`: Matriz (steps × n_features) de observables medidos
- `inputs`: Vector de inputs (0s y 1s)
- `window`: Ventana temporal para considerar historia

# Retorna
Matriz de separabilidad entre pares de tiempos
"""
function separation_kernel(trajectory::Matrix, inputs::Vector; window::Int=1)
    steps = size(trajectory, 1)
    kernel = zeros(steps, steps)
    
    for i in 1:steps
        for j in i:steps
            # Solo calcular si hay diferencia en los inputs en la ventana
            if i > window && j > window
                # Comparar secuencias de inputs
                seq_i = inputs[max(1, i-window):i]
                seq_j = inputs[max(1, j-window):j]
                
                if seq_i != seq_j
                    # Distancia entre estados del reservorio
                    kernel[i,j] = euclidean_distance(trajectory[i,:], trajectory[j,:])
                    kernel[j,i] = kernel[i,j]
                end
            end
        end
    end
    
    return kernel
end

# ==============================================================================
# CAPACIDAD DE SEPARACIÓN (Separation Capacity)
# ==============================================================================

"""
    separation_capacity(trajectory_A::Matrix, trajectory_B::Matrix, 
                       inputs::Vector; threshold::Float64=0.1)

Calcula la capacidad de separación entre dos trayectorias iniciadas 
desde estados diferentes con los mismos inputs.

# Argumentos
- `trajectory_A`: Trayectoria desde estado inicial A
- `trajectory_B`: Trayectoria desde estado inicial B  
- `inputs`: Secuencia de inputs compartida
- `threshold`: Umbral para considerar separación significativa

# Retorna
- Fracción de pasos donde la separación excede el umbral
"""
function separation_capacity(trajectory_A::Matrix, trajectory_B::Matrix, 
                             inputs::Vector; threshold::Float64=0.1)
    steps = size(trajectory_A, 1)
    significant_separations = 0
    
    for k in 2:steps
        if inputs[k] != inputs[k-1]  # Solo en transiciones
            dist = euclidean_distance(trajectory_A[k,:], trajectory_B[k,:])
            if dist > threshold
                significant_separations += 1
            end
        end
    end
    
    # Número total de transiciones
    n_transitions = sum(inputs[2:end] .!= inputs[1:end-1])
    
    return n_transitions > 0 ? significant_separations / n_transitions : 0.0
end

# ==============================================================================
# SEPARABILIDAD PROMEDIO EN VENTANAS
# ==============================================================================

"""
    windowed_separation(trajectory::Matrix, inputs::Vector, window_size::Int=10)

Calcula la separabilidad promedio en ventanas deslizantes.
Útil para visualizar evolución temporal de la separabilidad.

# Retorna
Vector con separabilidad promedio en cada ventana
"""
function windowed_separation(trajectory::Matrix, inputs::Vector, window_size::Int=10)
    steps = size(trajectory, 1)
    n_windows = steps - window_size + 1
    sep_values = zeros(n_windows)
    
    for i in 1:n_windows
        window_traj = trajectory[i:i+window_size-1, :]
        window_inp = inputs[i:i+window_size-1]
        
        # Calcular separación promedio en esta ventana
        sep_sum = 0.0
        count = 0
        for k in 2:window_size
            if window_inp[k] != window_inp[k-1]
                sep_sum += euclidean_distance(window_traj[k,:], window_traj[k-1,:])
                count += 1
            end
        end
        
        sep_values[i] = count > 0 ? sep_sum / count : 0.0
    end
    
    return sep_values
end

# ==============================================================================
# ÍNDICE DE SEPARABILIDAD GLOBAL
# ==============================================================================

"""
    global_separation_index(separation_distances::Vector)

Calcula un índice global de separabilidad a partir del vector de distancias.

# Retorna
Tupla (mean, std, max) de las separaciones
"""
function global_separation_index(separation_distances::Vector)
    non_zero = separation_distances[separation_distances .> 0]
    
    if length(non_zero) == 0
        return (0.0, 0.0, 0.0)
    end
    
    mean_sep = mean(non_zero)
    std_sep = std(non_zero)
    max_sep = maximum(non_zero)
    
    return (mean_sep, std_sep, max_sep)
end

# ==============================================================================
# FUNCIÓN AUXILIAR: CALCULAR TODAS LAS MÉTRICAS
# ==============================================================================

"""
    compute_all_separation_metrics(hist_A, hist_B, separation_dist, inputs)

Calcula todas las métricas de separabilidad disponibles.

# Retorna
Dict con todas las métricas calculadas
"""
function compute_all_separation_metrics(hist_A::Matrix, hist_B::Matrix, 
                                       separation_dist::Vector, inputs::Vector)
    
    # Índices globales
    mean_sep, std_sep, max_sep = global_separation_index(separation_dist)
    
    # Capacidad de separación
    sep_capacity = separation_capacity(hist_A, hist_B, inputs)
    
    # Kernel de separabilidad (solo para trayectoria A)
    sep_kernel = separation_kernel(hist_A, inputs, window=3)
    
    # Separabilidad por ventanas
    windowed_sep = windowed_separation(hist_A, inputs, 10)
    
    return Dict(
        "mean_separation" => mean_sep,
        "std_separation" => std_sep,
        "max_separation" => max_sep,
        "separation_capacity" => sep_capacity,
        "separation_kernel" => sep_kernel,
        "windowed_separation" => windowed_sep,
        "instant_separations" => separation_dist
    )
end