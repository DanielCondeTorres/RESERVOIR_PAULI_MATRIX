using JLD2

# La ruta exacta que falló en tu log
ruta = "/Users/danielcondetorres/Desktop/IBM_EXAMENES/TFM/RESERVOIR_PAULI_MATRIX/Code/Input_Data/6_1_2_all_zeros_12345.jld2"

println("📂 Abriendo archivo...")
data = load(ruta)

println("\n🔑 Claves en la RAÍZ del archivo:")
println(keys(data))

if haskey(data, "meta_data_dict")
    println("\n📋 Claves dentro de 'meta_data_dict':")
    # Imprimimos cada clave para ver si tiene otro nombre (ej: "s_k", "forcing", etc.)
    for k in keys(data["meta_data_dict"])
        println(" - $k")
    end
else
    println("\n❌ No existe 'meta_data_dict' en este archivo.")
end