using HDF5
using CSV
using DataFrames

# --- 1. CONFIGURACIÓN DE RUTAS ---
archivo_jld2 = expanduser("~/Desktop/IBM_EXAMENES/TFM/RESERVOIR_PAULI_MATRIX/Code/Input_Data/results_4_STM.jld2")
directorio_salida = expanduser("~/Desktop/IBM_EXAMENES/TFM/RESERVOIR_PAULI_MATRIX/Code/SEE_NATHAN/")
archivo_csv_final = joinpath(directorio_salida, "MASTER_DATA_CONSOLIDATED.csv")

# Crear la carpeta por si acaso
mkpath(directorio_salida)

println("🚀 Iniciando Extracción Maestra...")
println("📂 Origen: $archivo_jld2")

# Lista para acumular los datos
datos_acumulados = []

# Función para procesar y acumular
function acumular_numeros(datos, etiqueta)
    if datos isa AbstractArray{<:Number} || datos isa AbstractArray{ComplexF64}
        val_vec = vec(collect(datos))
        for (idx, v) in enumerate(val_vec)
            push!(datos_acumulados, (
                Ruta = etiqueta,
                Indice = idx,
                Real = real(v),
                Imag = imag(v)
            ))
        end
        return true
    end
    return false
end

# Función recursiva profunda
function explorar_y_consolidar(obj, nombre, h5file)
    # A. Si es un Dataset
    if obj isa HDF5.Dataset
        contenido = try read(obj) catch; return end
        
        # 1. ¿Son números directos?
        if acumular_numeros(contenido, nombre)
            return
        end

        # 2. ¿Es una referencia?
        if contenido isa HDF5.Reference
            try explorar_y_consolidar(h5file[contenido], nombre * "_ref", h5file) catch; end
        
        # 3. ¿Es una lista de referencias? (Diccionarios JLD2)
        elseif contenido isa AbstractArray{HDF5.Reference}
            for (i, r) in enumerate(contenido)
                try explorar_y_consolidar(h5file[r], nombre * "_i$i", h5file) catch; end
            end
            
        # 4. ¿Es un NamedTuple (como el kvvec)?
        elseif contenido isa NamedTuple
            for campo in keys(contenido)
                val = getfield(contenido, campo)
                if val isa HDF5.Reference
                    try explorar_y_consolidar(h5file[val], nombre * "_$campo", h5file) catch; end
                end
            end
        end

    # B. Si es un Grupo
    elseif obj isa HDF5.Group
        for k in keys(obj)
            if !startswith(k, "_")
                explorar_y_consolidar(obj[k], nombre * "/" * k, h5file)
            end
        end
    end
end

# --- 2. EJECUCIÓN ---
h5open(archivo_jld2, "r") do file
    for k in keys(file)
        if !startswith(k, "_")
            println("🔍 Escaneando rama: $k")
            explorar_y_consolidar(file[k], k, file)
        end
    end
end

println("\n------------------------------------------------")
total_puntos = length(datos_acumulados)
println("📊 Puntos de datos recolectados: $total_puntos")

if total_puntos > 0
    println("💾 Escribiendo archivo CSV en: $archivo_csv_final")
    df_final = DataFrame(datos_acumulados)
    CSV.write(archivo_csv_final, df_final)
    println("✅ ¡ÉXITO! Archivo creado correctamente.")
else
    println("❌ ERROR: No se encontró ningún dato numérico. El CSV no se creó.")
end