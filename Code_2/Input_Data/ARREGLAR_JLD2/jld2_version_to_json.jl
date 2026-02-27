import Pkg
Pkg.add("JSON")

using JLD2
using JSON

# Rutas de entrada y salida
ruta_jld2 = "/Users/danielcondetorres/Desktop/IBM_EXAMENES/TFM/RESERVOIR_PAULI_MATRIX/Code/Input_Data/results_clean.jld2"
ruta_json = "/Users/danielcondetorres/Desktop/IBM_EXAMENES/TFM/RESERVOIR_PAULI_MATRIX/Code/Input_Data/results.json"

println("🛡️ Cargando archivo JLD2...")

# Función para resolver referencias y valores simples
function resolver_valor(obj)
    if obj isa AbstractArray || obj isa Dict || obj isa Number || obj isa String || obj isa Bool
        return obj
    else
        # Convertir tuplas de enteros como claves a strings, porque JSON no admite tuplas como keys
        if obj isa Tuple
            return string(obj)
        else
            return string(obj)  # todo lo demás se convierte a string
        end
    end
end

# Abrimos JLD2 y construimos un diccionario limpio
datos_limpios = Dict()

jldopen(ruta_jld2, "r") do f
    for k in keys(f)
        try
            val = f[k]

            if val isa Dict
                # Convertimos todas las claves y valores a algo JSON-safe
                datos_limpios[k] = Dict(
                    resolver_valor(kk) => resolver_valor(vv) for (kk,vv) in val
                )
            else
                datos_limpios[k] = resolver_valor(val)
            end
        catch e
            println("⚠️ No se pudo leer la clave $k: $e")
        end
    end
end

println("💾 Guardando datos en JSON...")

open(ruta_json, "w") do io
    JSON.print(io, datos_limpios)
end

println("✅ Exportación completa a JSON: $ruta_json")

