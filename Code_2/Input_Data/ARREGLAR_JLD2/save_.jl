import Pkg
Pkg.add("JLD2")

using JLD2

# ---------------------------
# Rutas
# ---------------------------
ruta_original = "/Users/danielcondetorres/Desktop/IBM_EXAMENES/TFM/RESERVOIR_PAULI_MATRIX/Code/Input_Data/results_4_STM.jld2"
ruta_limpia   = "/Users/danielcondetorres/Desktop/IBM_EXAMENES/TFM/RESERVOIR_PAULI_MATRIX/Code/Input_Data/results_clean.jld2"

println("🛡️ Abriendo archivo original en Julia 1.8.5...")

# ---------------------------
# Abrir archivo original
# ---------------------------
data = Dict()
jldopen(ruta_original, "r") do f
    for k in keys(f)
        try
            # Intenta leer cada clave
            data[k] = f[k]
        catch e
            println("⚠️ No se pudo leer la clave $k, se omite.")
        end
    end
end

println("✅ Extracción completa. Claves disponibles: ", keys(data))

# ---------------------------
# Guardar archivo limpio
# ---------------------------
@save ruta_limpia data

println("🏁 Archivo limpio guardado en: $ruta_limpia")
println("Ahora puedes abrirlo en Julia 1.11 sin problemas de TypeVar.")

