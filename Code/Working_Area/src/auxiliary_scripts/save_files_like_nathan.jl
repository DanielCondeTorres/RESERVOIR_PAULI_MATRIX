using JLD2, FileIO

function save_qrc_results_jld2(N, steps, T_evol, h_val, inputs, expect_dict, folder_path)
    # 1. Convertimos la ruta a absoluta para evitar líos con el "../"
    # Esto asegura que se cree fuera de 'Working_Area' si usas "../"
    full_folder_path = abspath(joinpath(@__DIR__, "../../", folder_path))
    
    # --- LA CORRECCIÓN: 'isdir' en lugar de 'ispread' ---
    if !isdir(full_folder_path)
        mkpath(full_folder_path)
        println("📂 Carpeta creada en: $full_folder_path")
    end

    println("\n💾 Guardando datos en JLD2...")
    
    meta_data_dict = Dict(
        "n_qubits" => N,
        "n_steps"  => steps,
        "T_evol"   => T_evol,
        "h_val"    => h_val,
        "inputs"   => inputs,
        "method"   => "RK4_Matrix_Density_Full_Z_ZZ"
    )

    filename = "nathan_foundation_$(N)_steps$(steps)_full.jld2"
    full_file_path = joinpath(full_folder_path, filename)
    
    jldsave(full_file_path; meta_data_dict, expect_dict)
    println("✅ Archivo JLD2 guardado en: $full_file_path")
    
    return full_folder_path # Devolvemos la ruta para usarla en la gráfica
end