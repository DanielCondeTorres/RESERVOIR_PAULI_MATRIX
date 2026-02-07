"""
    separability_plots.jl

Módulo para visualizar métricas de separabilidad en QRC.
"""

using Plots
using Statistics

# ==============================================================================
# PLOT 1: EVOLUCIÓN TEMPORAL DE LA SEPARACIÓN
# ==============================================================================

"""
    plot_separation_evolution(separation_dist, inputs, save_path)

Gráfico de la evolución temporal de la distancia de separación.
Muestra cuándo el sistema separa estados tras cambios de input.
"""
function plot_separation_evolution(separation_dist::Vector, inputs::Vector, save_path::String)
    steps = length(separation_dist)
    
    # Identificar transiciones de input
    transitions = [0; diff(inputs)]
    transition_indices = findall(x -> x != 0, transitions)
    
    p = plot(1:steps, separation_dist, 
             label="Distancia de Separación",
             linewidth=2,
             color=:blue,
             xlabel="Paso Temporal",
             ylabel="Distancia (||Δx||₂)",
             title="Evolución de la Separabilidad",
             legend=:topright,
             size=(1000, 400),
             grid=true,
             gridalpha=0.3)
    
    # Marcar transiciones de input
    vline!(p, transition_indices, 
           label="Transiciones de Input",
           color=:red,
           linestyle=:dash,
           alpha=0.3)
    
    # Línea de referencia (separación promedio)
    non_zero = separation_dist[separation_dist .> 0]
    if length(non_zero) > 0
        hline!(p, [mean(non_zero)],
               label="Separación Promedio",
               color=:green,
               linestyle=:dot,
               linewidth=2)
    end
    
    savefig(p, joinpath(save_path, "separation_evolution.png"))
    return p
end

# ==============================================================================
# PLOT 2: HEATMAP DEL KERNEL DE SEPARABILIDAD
# ==============================================================================

"""
    plot_separation_kernel(kernel_matrix, save_path)

Heatmap del kernel de separabilidad SP(i,j).
Muestra qué pares de tiempos tienen estados más separados.
"""
function plot_separation_kernel(kernel_matrix::Matrix, save_path::String)
    p = heatmap(kernel_matrix,
                c=:viridis,
                xlabel="Paso Temporal j",
                ylabel="Paso Temporal i",
                title="Kernel de Separabilidad SP(i,j)",
                colorbar_title="Distancia",
                size=(700, 600),
                aspect_ratio=:equal)
    
    savefig(p, joinpath(save_path, "separation_kernel.png"))
    return p
end

# ==============================================================================
# PLOT 3: HISTOGRAMA DE SEPARACIONES
# ==============================================================================

"""
    plot_separation_histogram(separation_dist, save_path)

Histograma de las distancias de separación.
Muestra la distribución de cuán separables son los estados.
"""
function plot_separation_histogram(separation_dist::Vector, save_path::String)
    non_zero = separation_dist[separation_dist .> 0]
    
    if length(non_zero) == 0
        @warn "No hay separaciones no-cero para graficar"
        return nothing
    end
    
    p = histogram(non_zero,
                  bins=30,
                  normalize=:probability,
                  xlabel="Distancia de Separación",
                  ylabel="Frecuencia Normalizada",
                  title="Distribución de Separabilidad",
                  label="",
                  color=:steelblue,
                  size=(800, 500),
                  grid=true)
    
    # Añadir líneas de estadísticas
    vline!(p, [mean(non_zero)],
           label="Media",
           color=:red,
           linewidth=2,
           linestyle=:dash)
    
    vline!(p, [median(non_zero)],
           label="Mediana",
           color=:green,
           linewidth=2,
           linestyle=:dot)
    
    savefig(p, joinpath(save_path, "separation_histogram.png"))
    return p
end

# ==============================================================================
# PLOT 4: SEPARABILIDAD POR VENTANAS
# ==============================================================================


# ==============================================================================
# PLOT 5: COMPARACIÓN DE TRAYECTORIAS A vs B
# ==============================================================================

"""
    plot_trajectory_comparison(hist_A, hist_B, qubit_idx, inputs, save_path)

Compara las trayectorias de un qubit específico desde estados iniciales diferentes.
Útil para visualizar el efecto de la separabilidad.
"""
function plot_trajectory_comparison(hist_A::Matrix, hist_B::Matrix, 
                                   qubit_idx::Int, inputs::Vector, save_path::String)
    steps = size(hist_A, 1)
    
    p = plot(1:steps, hist_A[:, qubit_idx],
             label="Trayectoria A (|000...⟩)",
             linewidth=2,
             color=:blue,
             xlabel="Paso Temporal",
             ylabel="⟨Zᵢ⟩",
             title="Comparación de Trayectorias - Qubit $qubit_idx",
             size=(1200, 500),
             grid=true)
    
    plot!(p, 1:steps, hist_B[:, qubit_idx],
          label="Trayectoria B (|111...⟩)",
          linewidth=2,
          color=:red,
          linestyle=:dash)
    
    # Marcar inputs en el fondo
    for k in 1:steps
        if inputs[k] == 1
            vspan!(p, [k-0.5, k+0.5], 
                   color=:gray, 
                   alpha=0.1, 
                   label="")
        end
    end
    
    savefig(p, joinpath(save_path, "trajectory_comparison_q$qubit_idx.png"))
    return p
end

# ==============================================================================
# PLOT 6: PANEL RESUMEN DE MÉTRICAS
# ==============================================================================

"""
    plot_metrics_summary(metrics_dict, save_path)

Panel resumen con todas las métricas de separabilidad calculadas.
"""
function plot_metrics_summary(metrics_dict::Dict, save_path::String)
    # Crear texto con las métricas
    text_content = """
    MÉTRICAS DE SEPARABILIDAD
    ========================
    
    Separación Media: $(round(metrics_dict["mean_separation"], digits=4))
    Desv. Estándar:   $(round(metrics_dict["std_separation"], digits=4))
    Separación Máxima: $(round(metrics_dict["max_separation"], digits=4))
    
    Capacidad de Separación: $(round(metrics_dict["separation_capacity"]*100, digits=2))%
    
    Número de Separaciones Significativas: $(sum(metrics_dict["instant_separations"] .> 0))
    """
    
    p = plot(annotation=(0.5, 0.5, text(text_content, 12, :left, :black, "Courier")),
             framestyle=:none,
             size=(600, 400),
             xlims=(0,1),
             ylims=(0,1),
             title="Resumen de Métricas")
    
    savefig(p, joinpath(save_path, "metrics_summary.png"))
    return p
end

# ==============================================================================
# FUNCIÓN PRINCIPAL: GENERAR TODOS LOS PLOTS
# ==============================================================================

"""
    generate_all_separability_plots(hist_A, hist_B, separation_dist, 
                                   inputs, metrics_dict, save_path)

Genera todos los plots de separabilidad y los guarda.
"""
function generate_all_separability_plots(hist_A::Matrix, hist_B::Matrix,
                                        separation_dist::Vector, inputs::Vector,
                                        metrics_dict::Dict, save_path::String)
    
    println("📊 Generando plots de separabilidad...")
    
    # 1. Evolución temporal
    plot_separation_evolution(separation_dist, inputs, save_path)
    println("  ✓ Evolución temporal")
    
    # 2. Kernel de separabilidad
    plot_separation_kernel(metrics_dict["separation_kernel"], save_path)
    println("  ✓ Kernel de separabilidad")
    
    # 3. Histograma
    plot_separation_histogram(separation_dist, save_path)
    println("  ✓ Histograma de separaciones")
    
    # 4. Ventanas deslizantes
    plot_windowed_separation(metrics_dict["windowed_separation"], 10, save_path)
    println("  ✓ Separabilidad por ventanas")
    
    # 5. Comparación de trayectorias (primeros 3 qubits)
    N = size(hist_A, 2)
    for q in 1:min(3, N)
        plot_trajectory_comparison(hist_A, hist_B, q, inputs, save_path)
    end
    println("  ✓ Comparación de trayectorias")
    
    # 6. Resumen de métricas
    plot_metrics_summary(metrics_dict, save_path)
    println("  ✓ Panel de resumen")
    
    println("✅ Todos los plots generados en: $save_path")
end










using StatsPlots # Necesario para la función boxplot

"""
    plot_separability_boxplots(dict_A, inputs, N, save_dir)

Crea una cuadrícula de boxplots comparando la respuesta de cada qubit 
ante los dos posibles valores de entrada (0 y 1).
"""
function plot_separability_boxplots(dict_A, inputs, N, save_dir; skip_steps=20)
    println("📦 Generando Boxplots de separabilidad (saltando $skip_steps pasos)...")

    z_labels = sort([l for l in keys(dict_A) if count(c -> c == 'Z', l) == 1])
    rows = Int(ceil(N/3))
    p = plot(layout=(rows, 3), size=(1000, 350*rows), title="Separability (after step $skip_steps)")

    for (i, lbl) in enumerate(z_labels)
        # Filtramos para ignorar los primeros pasos (transitorios)
        valid_indices = (skip_steps+1):length(inputs)
        
        filtered_inputs = inputs[valid_indices]
        filtered_data = dict_A[lbl][valid_indices]

        data0 = filtered_data[filtered_inputs .== 0]
        data1 = filtered_data[filtered_inputs .== 1]
        
        x_labels = [fill("Input 0", length(data0)); fill("Input 1", length(data1))]
        y_values = [data0; data1]

        boxplot!(p[i], x_labels, y_values, 
                 title="Qubit $i", ylabel="<$lbl>", 
                 legend=false, color=[:blue :red], alpha=0.6)
    end

    savefig(p, joinpath(save_dir, "Separability_Boxplots_Filtered.png"))
end