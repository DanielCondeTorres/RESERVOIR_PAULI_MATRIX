using JSON
using Plots

function guardar_graficas_individuales(ruta_archivo::String, experiment_path::String)
    # 1. Verificar archivo
    if !isfile(ruta_archivo)
        error("No encuentro el archivo JSON: ", ruta_archivo)
    end

    # 2. Configurar directorios
    output_dir = joinpath(experiment_path, "resultados_a_comparar")
    if !ispath(output_dir)
        mkpath(output_dir)
        println("📁 Carpeta creada: ", output_dir)
    end

    # 3. Cargar datos
    datos = JSON.parsefile(ruta_archivo)
    root = datos["data"]
    llave_simulacion = collect(keys(root))[1] 
    contenido = root[llave_simulacion]
    meta = contenido["meta"]
    expect = contenido["expect_dict"]
    
    s_vec = Float64.(meta["s_vec"])
    pasos = 1:length(s_vec)

    println("🚀 Generando gráficas individuales en: ", output_dir)

    # 4. Iterar sobre cada observable
    for k in keys(expect)
        valores = Float64.(expect[k])
        
        # Opcional: saltar si es todo ceros para no generar basura
        if sum(abs.(valores)) < 1e-9
            continue
        end

        # Crear un plot de 2 filas para este observable específico
        p = plot(layout=(2,1), size=(800, 600), dpi=150, titlefont=10)

        # Gráfico superior: Entrada
        plot!(p[1], pasos, s_vec, 
              label="Input (s_vec)", color=:black, lw=2, linetype=:steppost,
              ylabel="Amplitud", title="Observable: $k (N=$(meta["N"]), g=$(meta["g"]))")

        # Gráfico inferior: El observable específico
        plot!(p[2], pasos, valores, 
              label=k, color=:red, lw=2.5, marker=:circle,
              xlabel="Steps", ylabel="⟨$k⟩", grid=true)
        ylims!(p, -1.1, 1.1)
        # 5. Guardar la imagen con el nombre del observable
        nombre_archivo = joinpath(output_dir, "$k.png")
        savefig(p, nombre_archivo)
        println("✅ Guardado: $k.png")
    end

    println("\n✨ ¡Proceso finalizado! Revisa la carpeta 'resultados_a_comparar'.")
end