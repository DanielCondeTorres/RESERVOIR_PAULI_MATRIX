function save_qrc_results_jld2(N, steps, T_evol, h_val, inputs, expect_dict, folder_path)
    # folder_path ahora usa la ruta directa que definiste arriba (Experiment_path)
    if !isdir(folder_path)
        mkpath(folder_path)
        println("📂 Carpeta creada en: $folder_path")
    end

    println("\n💾 Guardando datos en JLD2...")
    
    meta_data_dict = Dict(
        "n_qubits" => N,
        "n_steps"  => steps,
        "T_evol"   => T_evol,
        "h_val"    => h_val,
        "inputs"   => inputs,
        "method"   => "Pauli_Propagation_RK4_Full"
    )

    filename = "nathan_foundation_$(N)_steps$(steps)_full.jld2"
    full_file_path = joinpath(folder_path, filename) # Guardado directo y seguro
    
    jldsave(full_file_path; meta_data_dict, expect_dict)
    
    println("✅ ¡ARCHIVO GUARDADO FÍSICAMENTE AQUÍ! 👇")
    println("📍 $full_file_path")
    
    return folder_path
end