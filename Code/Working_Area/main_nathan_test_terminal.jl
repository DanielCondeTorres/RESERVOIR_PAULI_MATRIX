# ==============================================================================
# 1. IMPORTACIONES (¡ESTO ERA LO QUE FALTABA!)
# ==============================================================================
using LinearAlgebra
using Plots
using JLD2
using FileIO
using Random

# ==============================================================================
# 2. INCLUDES (Asegúrate de que las rutas son correctas en tu PC)
# ==============================================================================
# Ajusta la ruta "../src/..." según donde guardes este archivo.
# Si estás en Code/Working_Area, esto debería funcionar:
include("src/operator_terms/pauli.jl")      
include("src/utils/pauli_algebra.jl")       
include("src/utils/dynamics.jl")
include("src/utils/injection.jl")
include("src/utils/initial_state.jl") 
include("src/utils/quantum_channels.jl") 
include("src/utils/shot_noise.jl") 

# ==============================================================================
# 3. FUNCIONES AUXILIARES
# ==============================================================================


function parse_nathan_key(key_str::String)
    # Busca índices donde hay 'Z'
    active_indices = Int[]
    for (i, char) in enumerate(key_str)
        if char == 'Z'
            push!(active_indices, i - 1) # Nathan 1-based -> Tu código 0-based
        end
    end
    return active_indices
end

# ==============================================================================
# 4. FUNCIÓN PRINCIPAL
# ==============================================================================
function run_correlation_validation()
    println("=== VALIDACIÓN DE CORRELACIONES (Z1 Z2) ===")
    
    # --- A. CONFIGURACIÓN ---
    # Índices de los Qubits a correlacionar (0-based)
    # Ejemplo: Qubit 1 y Qubit 2 (índices 1 y 2)
    target_indices = [1, 2] 
    
    # Construimos la máscara combinada (Bitwise OR)
    # Esto crea un operador que es Z en pos 1 Y Z en pos 2 simultáneamente
    target_z_mask = 0
    for idx in target_indices
        target_z_mask |= (1 << idx)
    end
    
    # --- B. CARGAR DATOS ---
    filename = "6_3_2_all_zeros_12345.jld2" # Archivo con datos Z
    
    # Búsqueda robusta del archivo
    path = ""
    possible_paths = [
        filename,
        joinpath("Input_Data", filename),
        "/Users/danielcondetorres/Desktop/IBM_EXAMENES/TFM/RESERVOIR_PAULI_MATRIX/Code/Input_Data/" * filename,
        joinpath(@__DIR__, "Input_Data", filename)
    ]
    for p in possible_paths; if isfile(p); path = p; break; end; end
    
    if path == ""
        println("⚠️ No encuentro el archivo Z ($filename).")
        println("   Usaré el archivo X (6_1_2) si existe, pero NO podré comparar con Nathan.")
        path = replace(possible_paths[3], "6_3_2" => "6_1_2")
        if !isfile(path); error("❌ No encuentro ningún archivo de datos."); end
    end
    println("📂 Cargando: $path")
    
    data = load(path)
    meta = data["meta_data_dict"]
    expect_dict = data["expect_dict"]
    
    # Buscar si Nathan tiene esta correlación calculada
    found_key = ""
    nathan_series = nothing
    
    println("🔎 Buscando clave para índices $target_indices...")
    for (key, val) in expect_dict
        if key isa String
            indices_in_key = parse_nathan_key(key)
            # Comparamos conjuntos de índices (orden no importa)
            if Set(indices_in_key) == Set(target_indices)
                found_key = key
                nathan_series = real.(val)
                println("   ✅ ¡Encontrada! Clave: '$key'")
                break
            end
        end
    end
    
    if found_key == ""
        println("⚠️ Nathan no guardó esta correlación específica. Solo simularé tu parte.")
        nathan_series = zeros(Int(meta["num_steps"])) # Dummy
    end

    # --- C. SIMULACIÓN ---
    J = Float64.(meta["Jvec"]) 
    h = Float64(meta["h"])
    g = Float64(meta["g"])
    inputs_s = Float64.(meta["s_vec"])
    N = Int(meta["N"])
    
    # Física: Ising Estándar (ZZ + X) -> La que funciona para inyección Z
    H_evol = Operator()
    for i in 0:(N-2); H_evol[PauliString(0, (1<<i)|(1<<(i+1)))] = -J[i+1]; end
    for i in 0:(N-1); H_evol[PauliString(1<<i, 0)] = -h; end

    rho = initial_state_all_zeros(N)
    scale_factor = 2.0^N
    dt = 0.01
    my_history = zeros(Float64, length(inputs_s))

    println("🚀 Simulando correlación <Z1 Z2>...")

    for k in 1:length(inputs_s)
        s_k = inputs_s[k]
        
        # Inyección Z
        rho = inject_state(rho, 0, 1.0 - 2.0*s_k, rx=2.0*sqrt(s_k*(1.0-s_k)))
        
        # Evolución (1000 pasos para T=10)
        for _ in 1:1000
            rho = step_rk4(rho, H_evol, dt)
            if length(rho) > 2500; truncate_operator!(rho, 2500); end
        end
        
        # Dephasing (Preserva Z)
        rho = apply_global_dephasingZ(rho, g)
        
        # --- MEDICIÓN DE CORRELACIÓN ---
        # Buscamos el coeficiente del operador Pauli Z1 Z2
        coeff = real(get(rho, PauliString(0, target_z_mask), 0.0im))
        
        # Escala y Ruido
        my_history[k] = apply_shot_noise(coeff, 1.5*10^5) * scale_factor
        
        if k % 100 == 0; print("\r   Paso $k..."); end
    end

    # --- D. RESULTADOS ---
    println("\n✅ Terminado.")
    p = plot(title="Correlación Z1 Z2", xlabel="Step", ylabel="<Z1 Z2>")
    
    if found_key != ""
        plot!(p, nathan_series[1:100], label="Nathan ($found_key)", lw=3, color=:blue)
        # Error solo si tenemos contra qué comparar
        err = norm(nathan_series - my_history) / norm(nathan_series)
        println("📊 Error Relativo: $(round(err * 100, digits=2)) %")
    end
    
    plot!(p, my_history[1:100], label="Tú (Simulación)", color=:red, ls=:dash, lw=2)
    
    # Crear carpeta Outputs si no existe
    out_dir = joinpath(@__DIR__, "Outputs")
    if !isdir(out_dir); mkpath(out_dir); end
    
    savefig(p, joinpath(out_dir, "correlation_validation2222.png"))
    println("🖼️ Guardado en Outputs/correlation_validation.png")
end

run_correlation_validation()
