using JLD2
using CSV
using DataFrames

# --- Configuración de Rutas ---
# Usamos expanduser para que entienda el símbolo "~"
archivo_entrada = expanduser("~/Desktop/IBM_EXAMENES/TFM/RESERVOIR_PAULI_MATRIX/Code/To_Nahan/resultado_INJ_DEPH_Z.jld2")
#archivo_entrada = expanduser("~/Desktop/IBM_EXAMENES/TFM/RESERVOIR_PAULI_MATRIX/archivos_comparar_nathan/resultado_simulacion_Z.jld2")
# Carpeta de salida (sin nombre de archivo, porque lo generaremos automático)
#carpeta_salida  = expanduser("~/Desktop/IBM_EXAMENES/TFM/RESERVOIR_PAULI_MATRIX/Code/csv_transformation_dani_to_nathan/")
carpeta_salida  = expanduser("~/Desktop/IBM_EXAMENES/TFM/RESERVOIR_PAULI_MATRIX/Code/csv_transformation_2/")
# Crear carpeta si no existe
mkpath(carpeta_salida)

# --- Cargar JLD2 ---
println("📂 Abriendo: $archivo_entrada")
contenido = load(archivo_entrada)
variables = keys(contenido)

println("Variables encontradas: $variables")
println("------------------------------------------------")

# --- Bucle: Procesar cada variable encontrada ---
for var_nombre in variables
    datos = contenido[var_nombre]
    println("🔄 Procesando variable: '$var_nombre' ...")
    
    # 1. Convertir a DataFrame según el tipo de dato
    df_final = try
        if datos isa DataFrame
            datos
        elseif datos isa Dict
            # Si es un Diccionario, lo convertimos a tabla Key-Value
            println("   ↳ Detectado Diccionario. Convirtiendo a formato Clave-Valor.")
            DataFrame(Clave = collect(keys(datos)), Valor = collect(values(datos)))
        else
            # Intento genérico para matrices u otros
            println("   ↳ Tipo de dato: $(typeof(datos)). Intentando conversión automática.")
            DataFrame(datos, :auto)
        end
    catch e
        println("   ❌ Error convirtiendo '$var_nombre'. Saltando...")
        continue
    end

    # 2. Generar nombre de archivo único
    # Ejemplo: 6_1_2_all_zeros_12_expect_dict.csv
    nombre_base = splitext(basename(archivo_entrada))[1]
    archivo_csv = joinpath(carpeta_salida, "$(nombre_base)_$(var_nombre).csv")

    # 3. Guardar
    CSV.write(archivo_csv, df_final)
    println("   ✅ Guardado: $archivo_csv")
end

println("------------------------------------------------")
println("¡Proceso terminado!")