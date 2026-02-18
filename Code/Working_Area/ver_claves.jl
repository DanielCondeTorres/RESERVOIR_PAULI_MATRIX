using JLD2

ruta_limpia = "/Users/danielcondetorres/Desktop/IBM_EXAMENES/TFM/RESERVOIR_PAULI_MATRIX/Code/Input_Data/results_clean.jld2"

println("🛡️ EXTRACCIÓN COMPLETA DE TODAS LAS CLAVES Y DATOS:")

function imprimir_objeto(obj, indent=0)
    pref = " " ^ (indent*2)  # sangría para jerarquía
    if obj isa AbstractArray
        println(pref, "Array de tamaño ", size(obj), ": ", obj[1:min(5,end)], " ...")
    elseif obj isa Dict
        println(pref, "Diccionario con ", length(obj), " claves")
        for (k,v) in obj
            println(pref, " 🔹 Clave: ", k)
            imprimir_objeto(v, indent+1)
        end
    else
        println(pref, typeof(obj), ": ", obj)
    end
end

jldopen(ruta_limpia, "r") do f
    for k in keys(f)
        println("\n📂 Clave raíz: ", k)
        try
            val = f[k]
            imprimir_objeto(val, 1)
        catch e
            println("⚠️ No se pudo leer el contenido de $k: ", e)
        end
    end
end

println("\n🏁 EXTRACCIÓN FINALIZADA.")
