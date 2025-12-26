using JLD2, FileIO

# ==============================================================================
# HELPER: GESTIÓN DE RUTAS (Privado)
# ==============================================================================
function _get_file_path(filename::String)
    root_folder = joinpath(@__DIR__, "..", "..", "..")
    path = joinpath(root_folder, "Input_Data", filename)
    
    if !isfile(path)
        path = joinpath(pwd(), "Input_Data", filename)
        if !isfile(path)
            error("❌ No encuentro el archivo: $filename")
        end
    end
    return path
end

# ==============================================================================
# HELPER: VERIFICADOR DE CLAVES (FLEXIBLE)
# ==============================================================================
"""
    matches_configuration(key_str, tipo_letra, qubits_objetivo) -> Bool
"""
function matches_configuration(key_str::String, tipo_letra::String, qubits_objetivo::Vector{Int})
    letra_buscada = uppercase(tipo_letra)[1] # 'Z', 'X' o 'Y'
    
    # IMPORTANTE: Asumimos que el String tiene índice base 1 (Julia).
    # Si tus Qubits son 0-based (0, 1, 2...), tienes dos opciones:
    # 1. Pasar [2] para pedir el Qubit 1.
    # 2. Descomentar la línea 'i_ajustado' abajo para restar 1 automáticamente.
    
    for (i, char_en_clave) in enumerate(key_str)
        
        # --- AJUSTE DE ÍNDICE ---
        # Si tu archivo es "012345" y tú pides Qubit 1:
        # La 'i' de Julia es 2. Si tus qubits_objetivo son [1], no coincidirá.
        # Si descomentas esto, mapeamos String[1]->Qubit[0], String[2]->Qubit[1]...
        
        qubit_actual = i - 1  # <--- ACTIVADO PARA TU CASO (1Z1111 -> Z es Qubit 1)
        
        if qubit_actual in qubits_objetivo
            # CASO 1: Es un qubit que queremos medir. Debe ser la letra (Z).
            if char_en_clave != letra_buscada
                return false
            end
        else
            # CASO 2: NO es un qubit objetivo. Debe ser Identidad.
            # ACEPTAMOS 'I' (letras) O '1' (números de Nathan)
            if char_en_clave != 'I' && char_en_clave != '1'
                return false
            end
        end
    end
    
    return true
end

# ==============================================================================
# FUNCIÓN PRINCIPAL: EXTRAER DATOS (UNIFICADA)
# ==============================================================================
function extract_nathan_data(filename::String, pauli_str::String, target_qubits::Vector{Int})
    println("\n=== 🔍 BUSCANDO DATOS ($pauli_str en Qubits $target_qubits) ===")
    
    file_path = _get_file_path(filename)
    data = load(file_path)
    
    if !haskey(data, "expect_dict"); error("❌ Falta 'expect_dict'"); end
    expect_dict = data["expect_dict"]
    
    serie_encontrada = nothing
    clave_encontrada = ""
    
    # --- DEBUG: Imprimir las primeras 3 claves para ver qué formato tienen ---
    println("   🕵️  (Debug) Primeras claves en el archivo:")
    for (i, k) in enumerate(keys(expect_dict))
        if i <= 3; println("       - '$k'"); end
    end
    # -----------------------------------------------------------------------

    for k in keys(expect_dict)
        if k isa String
            if matches_configuration(k, pauli_str, target_qubits)
                serie_encontrada = real.(expect_dict[k])
                clave_encontrada = k
                break 
            end
        end
    end

    if isnothing(serie_encontrada)
        println("⚠️  No se encontraron datos para: $pauli_str en $target_qubits")
        println("    Consejo: Revisa si tus índices empiezan en 0 o 1.")
        return nothing
    else
        println("   ✅ ¡Encontrada! Clave original: '$clave_encontrada'")
        return serie_encontrada
    end
end

# ==============================================================================
# FUNCIÓN METADATOS (Sin cambios)
# ==============================================================================
function extract_metadata(filename::String)
    println("\n=== ℹ️ LEYENDO METADATOS ===")
    file_path = _get_file_path(filename)
    data = load(file_path)
    if !haskey(data, "meta_data_dict"); error("❌ Falta 'meta_data_dict'"); end
    meta = data["meta_data_dict"]

    info = Dict(
        "N_qubits"  => Int(meta["N"]),             
        "J"         => Float64.(meta["Jvec"]),     
        "h"         => Float64(meta["h"]),         
        "g"         => Float64(meta["g"]),        
        "dt"        => Float64(meta["dt"]),        
        "num_steps" => Int(meta["num_steps"]),     
        "s_vec"     => Float64.(meta["s_vec"])     
    )
    println("   ✅ Metadatos cargados (N=$(info["N_qubits"]))")
    return info
end