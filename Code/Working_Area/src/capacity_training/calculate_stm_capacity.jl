using Statistics, LinearAlgebra, Printf

"""
    calculate_stm_capacity(reservoir_data::Dict, original_inputs::Vector, 
                           max_delay::Int=10, washout::Int=50)

Calcula la Capacidad de Memoria a Corto Plazo (STM) del reservorio.
Realiza el post-procesamiento, alineación de retardos (tau), división train/test
y cálculo de capacidad para cada tau.

Argumentos:
- `reservoir_data`: Diccionario con las trayectorias de los observables (ej: dict_B).
- `original_inputs`: Vector con la señal de entrada original.
- `max_delay`: Tau máximo a evaluar.
- `washout`: Pasos iniciales a descartar.

Retorna:
- `capacities`: Vector con la capacidad C para cada tau (0 a max_delay).
- `total_stm`: La suma total de capacidades.
"""
function calculate_stm_capacity(reservoir_data::Dict{String, Vector{Float64}}, 
                                original_inputs::Vector, 
                                max_delay::Int=10, 
                                washout::Int=50)

    println("🧠 Iniciando cálculo de STM (0 a $max_delay)...")
    
    # 1. PREPARACIÓN DE DATOS
    # Convertir inputs a Float si no lo son
    targets_float = Float64.(original_inputs)
    steps = length(targets_float)
    
    # Convertir Diccionario -> Matriz X [Steps x Features]
    # Ordenamos las keys alfabéticamente para asegurar determinismo
    feature_keys = sort(collect(keys(reservoir_data)))
    n_features = length(feature_keys)
    
    # Pre-allocating matriz
    X_reservoir = zeros(Float64, steps, n_features)
    for (j, key) in enumerate(feature_keys)
        # Aseguramos que la longitud coincida (por si acaso)
        data_col = reservoir_data[key]
        len_col = min(length(data_col), steps)
        X_reservoir[1:len_col, j] = data_col[1:len_col]
    end

    capacities = Float64[]
    train_ratio = 0.8  # 80% Train, 20% Test

    # 2. BUCLE DE RETARDOS (Tau)
    for tau in 0:max_delay
        # A. ALINEACIÓN TEMPORAL
        if tau == 0
            y_target = targets_float
            X_feats  = X_reservoir
        else
            # Predecir input hace 'tau' pasos:
            # El target pierde los últimos 'tau' valores
            y_target = targets_float[1 : end-tau]
            # El reservorio pierde los primeros 'tau' valores (empieza más tarde)
            X_feats  = X_reservoir[1+tau : end, :]
        end
        
        # B. LIMPIEZA DE WASHOUT
        # Quitamos el transitorio ANTES de dividir en Train/Test
        if length(y_target) <= washout
            println("⚠️ Advertencia: Tau $tau deja sin datos tras washout.")
            push!(capacities, 0.0)
            continue
        end

        X_clean = X_feats[washout+1:end, :]
        y_clean = y_target[washout+1:end]
        
        # C. SPLIT TRAIN / TEST
        n_split = floor(Int, length(y_clean) * train_ratio)
        if n_split < 10 # Seguridad mínima
            continue 
        end

        X_train = X_clean[1:n_split, :]
        y_train = y_clean[1:n_split]
        
        X_test  = X_clean[n_split+1:end, :]
        y_test  = y_clean[n_split+1:end]
        
        # D. ENTRENAMIENTO (Ridge/Linear)
        # Pasamos washout=0 porque ya limpiamos los datos en X_clean
        weights, _ = train_reservoir(X_train, y_train, 0)
        
        # E. PREDICCIÓN (Test)
        # Añadimos bias manual para test (hcat de 1s) porque 'weights' incluye el bias
        X_test_bias = hcat(X_test, ones(size(X_test, 1)))
        y_pred_test = X_test_bias * weights
        
        # F. CÁLCULO DE C
        C = calculate_capacity(y_test, y_pred_test)
        push!(capacities, C)
        
        # Print opcional para ver progreso
        # @printf("   Tau %2d -> C = %.4f\n", tau, C)
    end
    
    total_stm = sum(capacities)
    println("🏆 Capacidad Total STM = $(round(total_stm, digits=4))")
    
    return capacities, total_stm
end