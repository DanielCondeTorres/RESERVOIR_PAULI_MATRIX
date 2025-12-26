# Ubicación: main_preload.jl (en Working_Area)
using LinearAlgebra
using Plots
using JLD2
using FileIO

# --- CARGAR DEPENDENCIAS ---
include("src/utils/qrc_training.jl")       
include("src/visualization/plot_stm_capacity.jl") 
include("src/utils/data_loader.jl") 

function run_data_driven_stm()
    println("=== STM CAPACITY (Datos de Nathan) ===")
    
    # -----------------------------------------------------------
    # 1. LOCALIZACIÓN AUTOMÁTICA DE DATOS
    # -----------------------------------------------------------
    # Buscamos en la carpeta hermana (../Input_Data) según tu captura
    path_sibling = joinpath(@__DIR__, "..", "Input_Data")
    # Por si acaso, buscamos también en carpeta hija
    path_child = joinpath(@__DIR__, "Input_Data")
    
    base_dir = ""
    if isdir(path_sibling)
        base_dir = path_sibling
    elseif isdir(path_child)
        base_dir = path_child
    else
        println("\n❌ ERROR CRÍTICO: No encuentro la carpeta 'Input_Data'")
        println("   Busqué en: $(abspath(path_sibling))")
        println("   Y en:      $(abspath(path_child))")
        error("Por favor, mueve la carpeta Input_Data al lugar correcto.")
    end
    
    println("✅ Carpeta de datos detectada en: $(abspath(base_dir))")

    # -----------------------------------------------------------
    # 2. CARGA Y FUSIÓN (X, Y, Z)
    # -----------------------------------------------------------
    # Archivos con seed 12345
    files_to_combine = [
        joinpath(base_dir, "6_1_2_all_zeros_12345.jld2"), # 1 = X
        joinpath(base_dir, "6_2_2_all_zeros_12345.jld2"), # 2 = Y
        joinpath(base_dir, "6_3_2_all_zeros_12345.jld2")  # 3 = Z
    ]
    
    println("\n1. Combinando experimentos...")
    features_matrix, target_signal = combine_experiments(files_to_combine)
    
    n_steps, n_features = size(features_matrix)
    println("   -> Matriz final: $n_steps pasos x $n_features features.")
    
    # -----------------------------------------------------------
    # 3. ENTRENAMIENTO (CAPACIDAD vs DELAY)
    # -----------------------------------------------------------
    println("\n2. Entrenando y evaluando memoria...")
    
    washout = 200
    max_delay = 10
    capacities = zeros(Float64, max_delay)
    
    for tau in 1:max_delay
        # Definir Target (Input pasado) y Features (Estado presente)
        target_vec = target_signal[1:end-tau]
        feat_mat = features_matrix[tau+1:end, :]
        
        # Entrenar
        weights, preds = train_reservoir(feat_mat, target_vec, washout)
        
        # Evaluar
        real_target = target_vec[washout+1:end]
        real_preds = preds[washout+1:end]
        C = calculate_capacity(real_target, real_preds)
        
        capacities[tau] = C
        println("   Tau $tau: C = $(round(C, digits=4))")
    end
    
    # -----------------------------------------------------------
    # 4. GRAFICAR
    # -----------------------------------------------------------
    println("\n3. Generando gráfico...")
    plot_stm_capacity(capacities, max_delay, 6)
end

run_data_driven_stm()