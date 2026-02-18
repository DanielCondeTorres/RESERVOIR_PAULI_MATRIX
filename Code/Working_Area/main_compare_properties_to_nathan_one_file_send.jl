#LOAD MODULES
using LinearAlgebra
using Plots
using Random
using JLD2
using FileIO
using Statistics
using JSON
using Plots
const SCRIPT_DIR = @__DIR__

# ==============================================================================
# 1. PARAMETERS
# ==============================================================================
# Load from files
#INPUT_FILE = "6_3_2_all_zeros_12345.jld2"
# Parameters to calculate the capacity
# --- CARGA DE INPUTS (Buscando s_vec) ---
#input_file_path = isfile(INPUT_FILE) ? INPUT_FILE : joinpath(SCRIPT_DIR, "../Input_Data/", INPUT_FILE) #"../Input_Data/" is the path to INPUT_FILE
#if !isfile(input_file_path); error("❌ No se encuentra el archivo de inputs: $input_file_path"); end
#println(" Cargando inputs desde: $input_file_path")
#data_loaded = load(input_file_path)
#meta = data_loaded["meta_data_dict"]
#nathan_results = data_loaded["expect_dict"]
# End load from files (can be commented if you want to use random inputs)


#NEW INPUT FILES USING JSON (for better readability), BECAUSE PREVIOUS PROBLEMS

json_raw = "/Users/danielcondetorres/Desktop/IBM_EXAMENES/TFM/RESERVOIR_PAULI_MATRIX/Code/Input_Data/results.json"


#Save files:
Experiment_name = "Experiment_Nathan_dephasing_in_X_matrix_form" #Where i save the files
Experiment_path = joinpath(SCRIPT_DIR, "../$Experiment_name")
projective_mode = false # Not using it
if !isdir(Experiment_path); mkpath(Experiment_path); end
# --- SETUP SISTEMA (Usando parámetros del archivo) ---
#T_evol = 10.0        
#Js= 1.
#N = 6#Int(meta["N"])
#h_val =10.0# Float64(meta["h"])
#gamma = 1.0 #Float64(meta["g"]) 
#J_vec_file = nothing  # Vector{Float64}(meta["Jvec"])
#Random.seed!(1234)
#steps =  Int(meta["num_steps"])
#steps = 100
#inputs = rand(steps)
#inputs = Vector{Float64}(meta["s_vec"]) 
#n_substeps =100
#dt = 10*Js# meta["dt"]

datos = JSON.parsefile(json_raw)

# 2. Navegar hasta la meta-data
# Obtenemos la primera llave (ej: "g=0.3_seed=1...")
llave_raiz = collect(keys(datos["data"]))[1]
meta = datos["data"][llave_raiz]["meta"]

# 3. ASIGNACIÓN AUTOMÁTICA DE TUS VARIABLES
N = Int(meta["N"])               # Ahora será 4 automáticamente
h_val = Float64(meta["h"])       # 10.0
gamma = Float64(meta["g"])       # 0.3
T_evol = Float64(meta["Δt"])     # 10.0
steps = Int(meta["num_steps"])   # 10
inputs = Float64.(meta["s_vec"]) # [1.0, 1.0, 1.0, 1.0, 0.0, ...]

# Para el J_vec_file (los acoplamientos):
# Como Nathan guardó un Jdict, extraemos solo los valores numéricos
J_dict = meta["Jdict"]
J_vec_file = Float64.(collect(values(J_dict))) # Esto te da [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]



GATE_WHERE_DEPHASING_IS_APPLIED="Z"
println("\n LOADED PARAMETERS:")
println("-"^40)
println("N (Qubits):      $N")
println("h (Transversal): $h_val")
println("Gamma (g):       $gamma")
println("Steps (Pasos):   $steps")
#println("dt (Paso temp):  $dt")
# Verificamos si J_vec_file existe antes de pedirle la longitud
longitud_j = isnothing(J_vec_file) ? "0 (None/Random)" : length(J_vec_file)
println("Longitud Jvec:   $longitud_j")
println("Longitud s_vec:  $(length(inputs))")


# ==============================================================================
# 1. FUNCTIONS (EXTERNAL TO THE MAIN)
# ==============================================================================
function include_rel(path...)
    include(joinpath(SCRIPT_DIR, path...))
end
using Random
struct PauliString
    x_mask::Int
    z_mask::Int
end

const Operator = Dict{PauliString, ComplexF64}

# Function that transform operators to density matrix (in order to recicle previous work)
function operator_to_dense_matrix(H_dict::Operator, n_qubits::Int)
    dim = 2^n_qubits
    H_mat = zeros(ComplexF64, dim, dim)
    I2, X, Z, Y = [1 0; 0 1], [0 1; 1 0], [1 0; 0 -1], [0 -im; im 0]
    for (pauli, coeff) in H_dict
        # Empezamos con la matriz de identidad de 1x1
        term_mat = ComplexF64[1.0;;] 
        for k in 0:(n_qubits - 1)
            x = (pauli.x_mask >> k) & 1
            z = (pauli.z_mask >> k) & 1
            gate = I2
            if x == 1 && z == 0; gate = X
            elseif x == 0 && z == 1; gate = Z
            elseif x == 1 && z == 1; gate = Y
            end
            # kron(Qubits_posteriores, Qubits_anteriores)
            term_mat = kron(gate, term_mat)
        end
        H_mat += coeff * term_mat
    end
    return H_mat
end
# --- FUNCIÓN 1: Matriz Densa -> Diccionario de Paulis ---
function dense_to_operator(rho_mat::Matrix{ComplexF64}, n_qubits::Int)
    dim = 2^n_qubits
    rho_op = Operator()
    
    # Hay 4^n_qubits combinaciones posibles de Paulis
    for i in 0:(4^n_qubits - 1)
        x_mask = 0
        z_mask = 0
        temp_i = i
        
        # Construimos el PauliString a partir del índice i (base 4)
        for k in 0:(n_qubits - 1)
            pauli_type = temp_i % 4 # 0:I, 1:X, 2:Z, 3:Y
            if pauli_type == 1; x_mask |= (1 << k)      # X
            elseif pauli_type == 2; z_mask |= (1 << k)   # Z
            elseif pauli_type == 3; x_mask |= (1 << k); z_mask |= (1 << k) # Y
            end
            temp_i ÷= 4
        end
        
        p = PauliString(x_mask, z_mask)
        
        # MATEMÁTICA: c_p = Tr(rho * P) / 2^N
        # Reutilizamos tu función para crear la matriz de este Pauli específico
        p_mat = operator_to_dense_matrix(Dict(p => 1.0+0im), n_qubits)
        coeff = tr(rho_mat * p_mat) / dim
        
        # Solo guardamos términos significativos para ahorrar memoria
        if abs(coeff) > 1e-12
            rho_op[p] = coeff
        end
    end
    return rho_op
end

#PLOTTING ECHO AND SEPARABILITY
function plot_and_save_validation_full(dict_A, dict_B,  N, steps, save_dir; label_A="Tray. A", label_B="Tray. B")
    println("📊 Generando reporte detallado ($label_A vs $label_B)...")
    all_labels = collect(keys(dict_A)) 
    # 1. Detectar qué base estamos usando (X, Y o Z)
    basis_char = 'Z' # Default
    for char in ['X', 'Y', 'Z']
        if any(k -> contains(k, string(char)), all_labels)
            basis_char = char
            break
        end
    end
    b_str = string(basis_char)
    # 2. Crear carpetas con nombres dinámicos
    path_indiv = joinpath(save_dir, "Individual_$(b_str)")
    path_pairs = joinpath(save_dir, "Correlations_$(b_str)$(b_str)")
    path_triads = joinpath(save_dir, "Correlations_$(b_str)$(b_str)$(b_str)") # <--- NUEVA CARPETA
    mkpath(path_indiv); mkpath(path_pairs); mkpath(path_triads)
    # 3. Filtrar etiquetas
    # Individuales (1 letra), Pares (2 letras), Tríadas (3 letras)
    z_labels = sort([l for l in all_labels if count(c -> c == basis_char, l) == 1])
    zz_labels = sort([l for l in all_labels if count(c -> c == basis_char, l) == 2])
    zzz_labels = sort([l for l in all_labels if count(c -> c == basis_char, l) == 3]) 
    # 4. Gráficas individuales
    for lbl in z_labels
        idx = findfirst(isequal(basis_char), lbl) 
        p = plot(dict_A[lbl], label=label_A, c=:blue, lw=1.5, title="ESP: Qubit $idx ($b_str)")
        plot!(p, dict_B[lbl], label=label_B, c=:red, ls=:dash, lw=1.5)
        ylims!(p, -1.1, 1.1)
        savefig(p, joinpath(path_indiv, "ESP_$(b_str)$idx.png"))
    end

    # 5. Gráficas para pares
    for lbl in zz_labels
        indices = [i for (i, c) in enumerate(lbl) if c == basis_char]
        if length(indices) >= 2
            tag = "$(b_str)$(indices[1])$(b_str)$(indices[2])"
            p = plot(dict_A[lbl], label=label_A, c=:blue, lw=1.5, title="ESP Pair: $tag")
            plot!(p, dict_B[lbl], label=label_B, c=:red, ls=:dash, lw=1.5)
            ylims!(p, -1.1, 1.1)
            savefig(p, joinpath(path_pairs, "ESP_$tag.png"))
        end
    end

    # 6. Gráficas para Tríadas (Xi Xj Xk)
    for lbl in zzz_labels
        indices = [i for (i, c) in enumerate(lbl) if c == basis_char]
        if length(indices) >= 3
            tag = "$(b_str)$(indices[1])$(b_str)$(indices[2])$(b_str)$(indices[3])"
            p = plot(dict_A[lbl], label=label_A, c=:blue, lw=1.0, title="ESP Triad: $tag")
            plot!(p, dict_B[lbl], label=label_B, c=:red, ls=:dash, lw=1.0)
            ylims!(p, -1.1, 1.1)
            savefig(p, joinpath(path_triads, "ESP_$tag.png"))
        end
    end

    # 7. Resumen combinado
    p1 = plot(title="Sample Comparison")
    if !isempty(z_labels)
        sample_lbl = z_labels[1]
        p1 = plot(dict_A[sample_lbl], label=label_A, c=:blue, title="Sample ($sample_lbl)")
        plot!(p1, dict_B[sample_lbl], label=label_B, c=:red, ls=:dash)
        ylims!(p1, -1.1, 1.1)
    end
end
# Capacity 


# SAVE FILES
function save_qrc_results_jld2(N, steps, T_evol, h_val, inputs, expect_dict, folder_path)
    # 1. Convertimos la ruta a absoluta para evitar líos con el "../"
    # Esto asegura que se cree fuera de 'Working_Area' si usas "../"
    full_folder_path = abspath(joinpath(@__DIR__, "../../../", folder_path))
    
    # --- LA CORRECCIÓN: 'isdir' en lugar de 'ispread' ---
    if !isdir(full_folder_path)
        mkpath(full_folder_path)
        println("📂 Carpeta creada en: $full_folder_path")
    end

    println("\n💾 Guardando datos en JLD2...")
    
    meta_data_dict = Dict(
        "n_qubits" => N,
        "n_steps"  => steps,
        "T_evol"   => T_evol,
        "h_val"    => h_val,
        "inputs"   => inputs,
        "method"   => "RK4_Matrix_Density_Full_Z_ZZ"
    )

    filename = "nathan_foundation_$(N)_steps$(steps)_full.jld2"
    full_file_path = joinpath(full_folder_path, filename)
    
    jldsave(full_file_path; meta_data_dict, expect_dict)
    println("✅ Archivo JLD2 guardado en: $full_file_path")
    
    return full_folder_path # Devolvemos la ruta para usarla en la gráfica
end




using JSON
using Plots

function guardar_graficas_individuales(ruta_archivo::String, experiment_path::String)
    # 1. Verificar archivo
    if !isfile(ruta_archivo)
        error("No encuentro el archivo JSON: ", ruta_archivo)
    end

    # 2. Configurar directorios
    output_dir = joinpath(experiment_path, "resultados_a_comparar")
    if !ispath(output_dir)
        mkpath(output_dir)
        println("📁 Carpeta creada: ", output_dir)
    end

    # 3. Cargar datos
    datos = JSON.parsefile(ruta_archivo)
    root = datos["data"]
    llave_simulacion = collect(keys(root))[1] 
    contenido = root[llave_simulacion]
    meta = contenido["meta"]
    expect = contenido["expect_dict"]
    
    s_vec = Float64.(meta["s_vec"])
    pasos = 1:length(s_vec)

    println("🚀 Generando gráficas individuales en: ", output_dir)

    # 4. Iterar sobre cada observable
    for k in keys(expect)
        valores = Float64.(expect[k])
        
        # Opcional: saltar si es todo ceros para no generar basura
        if sum(abs.(valores)) < 1e-9
            continue
        end

        # Crear un plot de 2 filas para este observable específico
        p = plot(layout=(2,1), size=(800, 600), dpi=150, titlefont=10)

        # Gráfico superior: Entrada
        plot!(p[1], pasos, s_vec, 
              label="Input (s_vec)", color=:black, lw=2, linetype=:steppost,
              ylabel="Amplitud", title="Observable: $k (N=$(meta["N"]), g=$(meta["g"]))")

        # Gráfico inferior: El observable específico
        plot!(p[2], pasos, valores, 
              label=k, color=:red, lw=2.5, marker=:circle,
              xlabel="Steps", ylabel="⟨$k⟩", grid=true)
        ylims!(p, -1.1, 1.1)
        # 5. Guardar la imagen con el nombre del observable
        nombre_archivo = joinpath(output_dir, "$k.png")
        savefig(p, nombre_archivo)
        println("✅ Guardado: $k.png")
    end

    println("\n✨ ¡Proceso finalizado! Revisa la carpeta 'resultados_a_comparar'.")
end
# ========================================================
# EJEMPLO DE USO (Pega tu JSON aquí abajo entre las comillas)
# ========================================================

# ==============================================================================
# 2. END AUXILIAR FUNCTIONS
# ==============================================================================





# ==============================================================================
# 3. KEY FUNCTIONS
# ==============================================================================
# Function to create our Hamiltonian
function build_nathan_all_to_all_XX(n_qubits::Int, h::Float64, J_vec::Union{Vector{Float64}, Nothing}=nothing)
    n_pairs = Int(n_qubits * (n_qubits - 1) / 2)
    H = Operator()
    
    # --- LÓGICA DE SELECCIÓN DE J ---
    # Definimos Js explícitamente para claridad
    Js = 1.0 
    
    actual_J = if isnothing(J_vec) || length(J_vec) != n_pairs
        if isnothing(J_vec)
            println("🎲 J_vec is 'nothing'. Generating random J_ij en [-1, 1] (Js=$Js) para N=$n_qubits...")
        else
            println("⚠️ WARNING: J_vec length mismatch. Expected $n_pairs, got $(length(J_vec)).")
            println("🎲 Generating new random J_ij values...")
        end
        # Genera valores uniformes en [-Js, Js]
        # (rand(n_pairs) .- 0.5) * 2.0 * Js  <-- Esto da el rango [-1, 1] si Js=1
        (rand(n_pairs) .- 0.5)*Js
    else
        J_vec
    end

    # 1. Campo Transversal: -(h/2) * Z_j
    # Según Mujal, h=10Js y se aplica como h/2
    for j in 0:(n_qubits - 1)
        p_z = PauliString(0, 1 << j) 
        H[p_z] = get(H, p_z, 0.0im) + h / 2.0
    end

    # 2. Interacción All-to-All: J_ij * X_i * X_j
    idx = 1 
    for i in 0:(n_qubits - 1)
        for j in (i + 1):(n_qubits - 1)
            J_ij = actual_J[idx]
            p_xx = PauliString((1 << i) | (1 << j), 0)
            H[p_xx] = get(H, p_xx, 0.0im) + J_ij
            idx += 1
        end
    end
    return H
end

# Function to create all initial states:
function initial_state_all_zeros(n_qubits::Int)
    rho = Operator()
    # Coeficiente c_k = Tr(rho P) / 2^N. 
    # Para <Z>=1, el coeficiente es 1/2^N.
    norm_coeff = (1.0 / 2.0)^n_qubits + 0.0im
    # --- CORRECCIÓN AQUÍ: PARÉNTESIS AÑADIDOS ---
    limit = (1 << n_qubits) - 1
    for z_mask_val in 0:limit
        p = PauliString(0, z_mask_val)
        rho[p] = norm_coeff
    end 
    return rho
end
function initial_state_all_ones(n_qubits::Int)
    rho = Operator()
    # El coeficiente base es 1/2^N
    base_coeff = (1.0 / 2.0)^n_qubits + 0.0im
    # Recorremos todas las combinaciones posibles de Z (2^N)
    limit = (1 << n_qubits) - 1
    for z_mask in 0:limit
        # Contamos cuántos Z hay en esta máscara (cuántos bits a 1)
        num_zs = count_ones(z_mask)
        sign = (num_zs % 2 == 0) ? 1.0 : -1.0
        p = PauliString(0, z_mask)
        rho[p] = sign * base_coeff
    end
    return rho
end
# INJECTION 
function inject_state_EraseWrite_matrix(rho::Matrix{ComplexF64}, qubit_idx::Int, 
    rz::Float64, rx::Float64=0.0, ry::Float64=0.0)
    dim = size(rho, 1)
    N_qubits = Int(log2(dim))
    # 1. Matrices básicas
    I2 = [1.0+0im 0.0; 0.0 1.0]
    Zero_Proj = [1.0+0im 0.0; 0.0 0.0] # |0><0|
    Flip_Op   = [0.0+0im 1.0; 0.0 0.0] # |0><1| (Pone 1 en 0)
    # 2. Construir Operadores de Kraus para Resetear qubit k a |0>
    # M0 = I...|0><0|...I
    # M1 = I...|0><1|...I
    # (Usamos qubit_idx + 1 porque Julia indexa desde 1)
    k = qubit_idx + 1
    M0 = (k==1) ? Zero_Proj : I2
    M1 = (k==1) ? Flip_Op   : I2
    for i in 2:N_qubits
        M0 = kron(M0, (i==k) ? Zero_Proj : I2)
        M1 = kron(M1, (i==k) ? Flip_Op   : I2)
    end
    # Aplicamos el Reset: rho -> |0>_k<0| (x) Tr_k(rho)
    rho_reset = M0 * rho * M0' + M1 * rho * M1'
    # 3. Rotar al estado objetivo (rx, ry, rz)
    # Calculamos la unitaria que lleva |0> -> (rx, ry, rz)
    norm = sqrt(rx^2 + ry^2 + rz^2)
    if norm < 1e-9; 
        return rho_reset
    end # Si vector nulo, dejamos en |0>
    # Ángulos esféricos
    nx, ny, nz = rx/norm, ry/norm, rz/norm
    theta = acos(nz)       # Polar
    phi = atan(ny, nx)     # Azimutal
    # U_rot = Rz(phi) * Ry(theta)
    Rz = [exp(-im*phi/2) 0; 0 exp(im*phi/2)]
    Ry = [cos(theta/2) -sin(theta/2); sin(theta/2) cos(theta/2)]
    U_qubit = Rz * Ry
    # Expandimos a todo el sistema
    U_tot = (k==1) ? U_qubit : I2
    for i in 2:N_qubits
        U_tot = kron(U_tot, (i==k) ? U_qubit : I2)
    end
    return U_tot * rho_reset * U_tot'
end


function apply_global_dephasing(rho::Operator, g::Float64, axis::String)
    new_rho = Operator()
    decay_base = exp(-(g^2) / 2.0)
    axis_norm = uppercase(axis)
    
    # --- VARIABLES DE DEBUG ---
    total_terms = 0
    killed_terms = 0
    # --------------------------

    for (p, coeff) in rho
        if abs(coeff) < 1e-15; continue; end
        total_terms += 1
        
        n_anti = 0
        
        if axis_norm == "Z"
            # Si medimos Z, mueren los que tienen X (X e Y)
            n_anti = count_ones(p.x_mask)
            
        elseif axis_norm == "X"
            # Si medimos X, mueren los que tienen Z (Z e Y)
            n_anti = count_ones(p.z_mask)
            
        elseif axis_norm == "Y"
            n_anti = count_ones(p.x_mask ⊻ p.z_mask)
        end
        
        if n_anti > 0
            new_rho[p] = coeff * (decay_base ^ n_anti)
            killed_terms += 1
        else
            new_rho[p] = coeff
        end
    end
    # --- IMPRIMIR DEBUG ---
    # Esto imprimirá en cada paso qué está pasando.
    # Si ves que los números cambian al cambiar axis="X" a "Z", ¡funciona!
    if total_terms > 0
        println("🔎 DEBUG Dephasing [$axis_norm]: De $total_terms términos, se amortiguaron $killed_terms")
    end
    # ----------------------
    
    return new_rho
end




function build_full_pauli_matrix(pauli_matrix_where_measure::Matrix{ComplexF64}, k::Int, n_qubits::Int)
    # 1. Definimos la identidad local de 2x2
    Identiy = ComplexF64[1 0; 0 1]
    # Si k es 1, empezamos con pauli_matrix_where_measure; si no, con la Identidad.
    op_total = (k == 1) ? pauli_matrix_where_measure : Identiy
    # OJOOOOOOOOOO: K==1 Y NO K=0 POR SER JULIA, ADEMAS FUERA DEL BUCLE PARA EMPEZAR A CREAR EL OPERADOR!!!
    # 3. producto de Kronecker (⊗)
    for i in 2:n_qubits
        # Si el índice actual i coincide con el qubit objetivo k, aplicamos LA MATRIZ DE PAULI Q MEDIMOS
        # En caso contrario, aplicamos la identidad I!!
        next_gate = (i == k) ? pauli_matrix_where_measure : Identiy
        op_total = kron(op_total, next_gate)
    end
    return op_total
end
function apply_dephasing_matrix(rho::Matrix{ComplexF64}, g::Float64, axis::String)
    if g <= 1e-9; return rho; end
    n_qubits = Int(log2(size(rho, 1))) # or charge them in the function, but should be the same¿?
    # Nathan notes (1 - epsilon) = exp(-g^2/2)
    epsilon =1.0-exp(-(g^2) / 2.0)
    # Equation (13)
    term_a = 1.0 - (epsilon / 2.0)
    term_b  = epsilon / 2.0
    # Definir el operador local W (X, Y o Z)
    pauli_matrix_where_measure = if axis == "X" 
        ComplexF64[0 1; 1 0]
    elseif axis == "Y" 
        ComplexF64[0 -im; im 0]
    else 
        ComplexF64[1 0; 0 -1] 
    end
    rho_curr = copy(rho)

    # Aplicar el mapa para cada qubit k (Q recorre todos los 1-local W)
    for k in 1:n_qubits
        # Q_k es el operador W aplicado solo al qubit k
        Q_k = build_full_pauli_matrix(pauli_matrix_where_measure, k, n_qubits)
        # APLICACIÓN DE LA EQ (13): (1 - epsilon/2)*rho + (epsilon/2)*Q*rho*Q
        rho_curr = term_a * rho_curr + term_b * (Q_k * rho_curr * Q_k)
    end
    return rho_curr
end
# ==============================================================================
# 3. KEY FUNCTIONS
# ==============================================================================

# ==============================================================================
# 4. SIMULATION
# ==============================================================================
function run_task()
    println("Starting experiment: $Experiment_name")
    # Setup Hamiltonian and initial states
    H_op = build_nathan_all_to_all_XX(N, h_val,J_vec_file)
    H_dense = operator_to_dense_matrix(H_op, N)
    U_evol = exp(-im * H_dense * T_evol)
    U_evol_adj = adjoint(U_evol)
    rho_A = operator_to_dense_matrix(initial_state_all_zeros(N), N) # States all zeros
    rho_B = operator_to_dense_matrix(initial_state_all_ones(N), N) # States all ones
    rho_C = operator_to_dense_matrix(initial_state_all_zeros(N), N) # State all zeros to check the separability
    
    #Random.seed!(1234)
    #inputs = rand(steps)
    
    # In initial state C the first 30 steps are going to be the same as A, them random
    inputs_C = copy(inputs)
    step_perturbation = 30
    if steps > step_perturbation
        inputs_C[step_perturbation+1:end] = rand(steps - step_perturbation)
    end
    ############################
    # OBSERVABLES TO MEASURE....
    ############################
    #Observables X, Y, Z 
    obs_matrices = Dict{String, Matrix{ComplexF64}}()
    all_labels = String[]
    sx=[0. 1.; 1. 0.]; sy=[0. -im; im 0.]; sz=[1. 0.; 0. -1.]; id=[1. 0.; 0. 1.]
    for (b_char, b_op) in zip(['X','Y','Z'], [sx, sy, sz])
        b_str = string(b_char)
        # 1-body
        for i in 1:N
            lbl = ["1" for _ in 1:N]; lbl[i] = b_str; s_lbl = join(lbl)
            push!(all_labels, s_lbl)
            op = (i==1) ? b_op : id; for k in 2:N; op = kron(op, (k==i) ? b_op : id); end
            obs_matrices[s_lbl] = op
        end
        # 2-body
        for i in 1:N, j in (i+1):N
            lbl = ["1" for _ in 1:N]; lbl[i]=b_str; lbl[j]=b_str; s_lbl = join(lbl)
            push!(all_labels, s_lbl)
            op = (i==1) ? b_op : id; for k in 2:N; op = kron(op, (k==i||k==j) ? b_op : id); end
            obs_matrices[s_lbl] = op
        end
        # 3-body (TRÍADAS) 
        for i in 1:N, j in (i+1):N, k in (j+1):N
            lbl = ["1" for _ in 1:N]
            lbl[i] = b_str; lbl[j] = b_str; lbl[k] = b_str
            s_lbl = join(lbl)
            push!(all_labels, s_lbl)
            # Construcción del operador: Kronecker product
            op = (1==i || 1==j || 1==k) ? b_op : id
            for m in 2:N
                op = kron(op, (m==i || m==j || m==k) ? b_op : id)
            end
            obs_matrices[s_lbl] = op
        end
    end
    ############################
    # OBSERVABLES TO MEASURE....
    ############################
    # Initialization of our dictionary to save our results
    dict_A = Dict(lbl => zeros(Float64, steps) for lbl in all_labels)
    dict_B = Dict(lbl => zeros(Float64, steps) for lbl in all_labels)
    dict_C = Dict(lbl => zeros(Float64, steps) for lbl in all_labels)  
    # C. Dinamics evolution (the core of the code)
    # For each step
    for k in 1:steps
        s_k = Float64(inputs[k]) # we create our input at random
        rz = 1.0 - 2.0 * s_k # this is the way we put it into the bloch sphere
        # Use more ore less the same to C
        s_k_C = Float64(inputs_C[k])
        rz_C = 1.0 - 2.0 * s_k_C
        rx = 0.0; ry = 0.0 
        # In the first qubit wa clean the previous data  and inject sk
        rho_A = inject_state_EraseWrite_matrix(rho_A, 0, rz, rx, ry)
        rho_B = inject_state_EraseWrite_matrix(rho_B, 0, rz, rx, ry)
        rho_C = inject_state_EraseWrite_matrix(rho_C, 0, rz_C, rx, ry)
        #for _ in 1:n_substeps
            rho_A = U_evol * rho_A * U_evol_adj
            rho_B = U_evol * rho_B * U_evol_adj
            rho_C = U_evol * rho_C * U_evol_adj
        #end
        ###### DEPHASING WITH CONVERSIONS
        # DEPHASING
        #rho_A_op = dense_to_operator(rho_A, N)
        #rho_B_op = dense_to_operator(rho_B, N)
        #rho_C_op = dense_to_operator(rho_C, N)
        # 3. Ruido: Aquí ocurre la magia de Daniel/Mujal
        #rho_A_op_noisy = apply_global_dephasing(rho_A_op, gamma, GATE_WHERE_DEPHASING_IS_APPLIED)
        #rho_B_op_noisy = apply_global_dephasing(rho_B_op, gamma, GATE_WHERE_DEPHASING_IS_APPLIED)
        #rho_C_op_noisy = apply_global_dephasing(rho_C_op, gamma, GATE_WHERE_DEPHASING_IS_APPLIED)
        # 4. Regreso: De Paulis a Matriz (Tu función original)
        #rho_A = operator_to_dense_matrix(rho_A_op_noisy, N)
        #rho_B = operator_to_dense_matrix(rho_B_op_noisy, N)
        #rho_C = operator_to_dense_matrix(rho_C_op_noisy, N)
        # DEPHASING END WITH CONVERSIONS

        #Dephasing direct operations
        rho_A=apply_dephasing_matrix(rho_A, gamma, GATE_WHERE_DEPHASING_IS_APPLIED)
        rho_B=apply_dephasing_matrix(rho_B, gamma, GATE_WHERE_DEPHASING_IS_APPLIED)
        rho_C=apply_dephasing_matrix(rho_C, gamma, GATE_WHERE_DEPHASING_IS_APPLIED)
        # end direct dephasing operations
        for (lbl, op) in obs_matrices
            # MEASURE!!!!!
            val_A = real(tr(rho_A * op))
            val_B = real(tr(rho_B * op))
            val_C = real(tr(rho_C * op))
            # PUT VALUES IN THE DICTIONARY
            dict_A[lbl][k] = val_A
            dict_B[lbl][k] = val_B
            dict_C[lbl][k] = val_C
        end
    end
    # D. PLOTTING    
    for basis in ["X", "Y", "Z"]
        println("   🎨 Procesando Base $basis...")
        basis_dir = joinpath(Experiment_path, "Basis_$basis")
        if !isdir(basis_dir); mkpath(basis_dir); end
        sub_dict_A = filter(p -> contains(p.first, basis), dict_A)
        sub_dict_B = filter(p -> contains(p.first, basis), dict_B)
        sub_dict_C = filter(p -> contains(p.first, basis), dict_C)
        # CHECK ECHO PROPERTY
        try
            dir_AB = joinpath(basis_dir, "Check_ECHO_A_B")
            plot_and_save_validation_full(
                sub_dict_A, sub_dict_B,  N, steps, dir_AB;
                label_A="Tray. A", 
                label_B="Tray. B"
            )
            # 2. CHECK SEPARABILITY
            dir_AC = joinpath(basis_dir, "Check_separability_A_C")
            mkpath(dir_AC)            
            plot_and_save_validation_full(
                sub_dict_A, sub_dict_C, N, steps, dir_AC;
                label_A="Tray. A (Origintal)", 
                label_B="Tray. C (Perturbation)"
            )
        # Observe an error  
        catch e
            println("Error en plots $basis: $e")
            Base.showerror(stdout, e, catch_backtrace())
        end
    end    
    save_qrc_results_jld2(N, steps, T_evol, h_val, inputs, dict_A, Experiment_name)

end
# ¡¡LLAMA A LA FUNCIÓN!!
plt =  guardar_graficas_individuales(json_raw,Experiment_path)
display(plt)
run_task()