using Plots, Statistics, LinearAlgebra, Printf

# --- TUS FUNCIONES DE ENTRENAMIENTO ---
function train_reservoir(observables::Matrix{Float64}, targets::Vector{Float64}, washout::Int=20)
    # 1. Descartar Washout
    X = observables[washout+1:end, :]
    y = targets[washout+1:end]
    
    # 2. Añadir Bias (columna de 1s)
    rows = size(X, 1)
    X_bias = hcat(X, ones(rows))
    
    # 3. Resolver Mínimos Cuadrados: w = (X'X)^-1 X'y
    # Nota: Si X es singular, añadir una pequeña regularización (Ridge) ayuda: 
    # weights = (X_bias' * X_bias + 1e-6 * I) \ (X_bias' * y)
    weights = X_bias \ y
    
    # 4. Predicciones
    predictions = X_bias * weights
    return weights, predictions
end

function calculate_capacity(target::Vector{Float64}, prediction::Vector{Float64})
    n = min(length(target), length(prediction))
    y = target[end-n+1:end]
    y_pred = prediction[end-n+1:end]
    
    cv = cov(y, y_pred)
    v_y = var(y)
    v_pred = var(y_pred)
    
    if v_y < 1e-12 || v_pred < 1e-12; return 0.0; end
    
    return (cv^2) / (v_y * v_pred)
end

# --- NUEVA FUNCIÓN PARA GRAFICAR CAPACIDAD ---
function plot_memory_capacity(capacities::Vector{Float64}, max_delay::Int, total_stm::Float64, save_dir::String)
    taus = 0:max_delay
    
    # Asegurarnos de que las longitudes coincidan (a veces el loop empieza en 1 o en 0)
    if length(capacities) < length(taus)
        # Si falta tau=0, ajustamos
        taus = 1:max_delay 
    end

    p = bar(taus, capacities, 
        label="Memory Capacity", 
        xlabel="Delay (τ)", 
        ylabel="Capacity C(τ)",
        title="STM Capacity (Total = $(round(total_stm, digits=2)))",
        legend=:topright,
        color=:skyblue,
        ylims=(0, 1.1),
        lw=0, alpha=0.8
    )
    
    # Línea de tendencia suave
    plot!(p, taus, capacities, color=:blue, label="", lw=2)
    
    savefig(p, joinpath(save_dir, "STM_Capacity_Plot.png"))
    println("📊 Gráfica de capacidad guardada en: $save_dir")
end