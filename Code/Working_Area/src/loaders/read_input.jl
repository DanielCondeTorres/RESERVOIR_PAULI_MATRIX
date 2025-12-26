using JLD2
using FileIO

# Ajusta la ruta si es necesario para apuntar a uno de los archivos que sí tienes
file_path = joinpath(@__DIR__, "../../..", "Input_Data", "6_1_2_all_zeros_12345.jld2")

if isfile(file_path)
    println("📂 Abriendo archivo: $file_path")
    data = load(file_path)
    
    # 1. Ver claves principales
    println("\n--- CLAVES PRINCIPALES ---")
    println(keys(data))
    
    # 2. Ver claves dentro de meta_data_dict
    if haskey(data, "meta_data_dict")
        meta = data["meta_data_dict"]
        println("\n--- CLAVES EN META_DATA ---")
        println(keys(meta))
        
        # Intentar adivinar cuál es el input
        for k in keys(meta)
            val = meta[k]
            if val isa Vector && length(val) > 100
                println("   -> Posible candidato para 's': '$k' (Tipo: $(typeof(val)), Longitud: $(length(val)))")
            end
        end
    else
        println("\n⚠️ No se encontró 'meta_data_dict'")
    end
else
    println("❌ No encuentro el archivo para inspeccionar. Revisa la ruta.")
end