#LOAD MODULES
using LinearAlgebra
using Plots
using Random
using JLD2
using FileIO
using Statistics

# ==============================================================================
# 1. FUNCTIONS (EXTERNAL TO THE MAIN)
# ==============================================================================
const SCRIPT_DIR = @__DIR__
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

# Function to create our Hamiltonian
function build_nathan_all_to_all_XX(n_qubits::Int, h::Float64; seed::Int=12345)
    Random.seed!(seed)
    H = Operator()
    # 1. Término de Campo Transversal: -(h/2) * Z_j en cada sitio
    # En tu máscara: PauliString(x_mask, z_mask)
    for j in 0:(n_qubits - 1)
        p_z = PauliString(0, 1 << j) # Z en el sitio j
        H[p_z] = get(H, p_z, 0.0im) - (h)
    end
    # 2. Término de Interacción All-to-All: J_ij * X_i * X_j
    # Conectamos todos los pares posibles (i < j)
    for i in 0:(n_qubits - 1)
        for j in (i + 1):(n_qubits - 1)
            # J_ij aleatorio entre -0.5 y 0.5 (como indica Nathan en sus notas)
            J_ij = (rand() - 0.5)   #*4
            # X_i * X_j -> x_mask activo en bits i y j
            p_xx = PauliString((1 << i) | (1 << j), 0)
            # Sumamos al operador (por si acaso ya existiera el término)
            H[p_xx] = get(H, p_xx, 0.0im) + J_ij
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

function train_reservoir(observables::Matrix{Float64}, targets::Vector{Float64}, washout::Int=20)
    # 1. Descartar Washout
    X = observables[washout+1:end, :]
    y = targets[washout+1:end]
    # 2. Añadir Bias (columna de 1s)
    rows = size(X, 1)
    X_bias = hcat(X, ones(rows))
    # 3. Resolver Mínimos Cuadrados: w = (X'X)^-1 X'y
    # Nota: Si X es singular, añadir una pequeña regularización (Ridge) ayuda: 
    # weights = (X_bias' * X_bias + 1e-6 * I) \ (X_bias' * y)
    weights = X_bias \ y
    # 4. Predicciones
    predictions = X_bias * weights
    return weights, predictions
end

function calculate_capacity(target::Vector{Float64}, prediction::Vector{Float64})
    n = min(length(target), length(prediction))
    y = target[end-n+1:end]
    y_pred = prediction[end-n+1:end]
    cv = cov(y, y_pred)
    v_y = var(y)
    v_pred = var(y_pred)
    return (v_y < 1e-12 || v_pred < 1e-12) ? 0.0 : (cv^2) / (v_y * v_pred)
end


function calculate_stm_capacity(
    reservoir_data::Dict{String, Vector{Float64}}, 
    original_inputs::Vector, 
    max_delay::Int=15, 
    washout::Int=50
)

    targets_float = Float64.(original_inputs)
    steps = length(targets_float)

    # Convertir dict → matriz exactamente como la segunda versión
    feature_keys = sort(collect(keys(reservoir_data)))
    X_reservoir = hcat([reservoir_data[k] for k in feature_keys]...)

    capacities = Float64[]

    for tau in 0:max_delay

        # Alineación temporal
        if tau == 0
            y_target = targets_float
            X_feats  = X_reservoir
        else
            y_target = targets_float[1:end-tau]
            X_feats  = X_reservoir[1+tau:end, :]
        end

        # Si no hay suficientes datos tras washout
        if length(y_target) <= washout
            push!(capacities, 0.0)
            continue
        end

        # Limpieza washout
        X_clean = X_feats[washout+1:end, :]
        y_clean = y_target[washout+1:end]

        # Split 80/20
        n_split = floor(Int, length(y_clean) * 0.8)
        if n_split < 10
            push!(capacities, 0.0)
            continue
        end

        X_train = X_clean[1:n_split, :]
        y_train = y_clean[1:n_split]

        X_test  = X_clean[n_split+1:end, :]
        y_test  = y_clean[n_split+1:end]

        # Entrenamiento (washout=0 porque ya limpiamos)
        weights, _ = train_reservoir(X_train, y_train, 0)

        # Predicción test (añadiendo bias manual)
        X_test_bias = hcat(X_test, ones(size(X_test, 1)))
        y_pred_test = X_test_bias * weights

        push!(capacities, calculate_capacity(y_test, y_pred_test))
    end

    return capacities, sum(capacities)
end


function plot_memory_capacity(
    capacities::Vector{Float64}, 
    max_delay::Int, 
    total_stm::Float64, 
    save_dir::String
)

    # Igual que la segunda versión
    taus = 0:length(capacities)-1

    p = bar(
        taus,
        capacities,
        label="Memory Capacity",
        xlabel="Delay (τ)",
        ylabel="Capacity C(τ)",
        title="STM Capacity (Total = $(round(total_stm, digits=2)))",
        color=:skyblue,
        ylims=(0, 1.1),
        alpha=0.8
    )

    plot!(p, taus, capacities, color=:blue, label="", lw=2)

    savefig(p, joinpath(save_dir, "STM_Capacity_Plot.png"))
end


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

# ==============================================================================
# 2. END AUXILIAR FUNCTIONS
# ==============================================================================

# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================



# ==============================================================================
# 2. PARAMETERS
# ==============================================================================
N = 6
steps = 100
T_evol = 1.0        
h_val = 1.0
gamma = 0.04 # I am not using it right now      
n_substeps = 100
dt = T_evol / n_substeps
Experiment_name = "Experiment" #Where i save the files
Experiment_path = joinpath(SCRIPT_DIR, "../$Experiment_name")
projective_mode = false # Not using it
if !isdir(Experiment_path); mkpath(Experiment_path); end
# Parameters to calculate the capacity
WASHOUT = 0      
MAX_DELAY = 10     
TRAIN_RATIO = 0.9  

# ==============================================================================
# 4. SIMULATION
# ==============================================================================
function run_task()
    println("Starting experiment: $Experiment_name")
    # Setup Hamiltonian and initial states
    H_op = build_nathan_all_to_all_XX(N, h_val)
    H_dense = operator_to_dense_matrix(H_op, N)
    U_evol = exp(-im * H_dense * T_evol)
    U_evol_adj = adjoint(U_evol)
    rho_A = operator_to_dense_matrix(initial_state_all_zeros(N), N) # States all zeros
    rho_B = operator_to_dense_matrix(initial_state_all_ones(N), N) # States all ones
    rho_C = operator_to_dense_matrix(initial_state_all_zeros(N), N) # State all zeros to check the separability
    
    Random.seed!(1234)
    inputs = rand(steps)
    
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
        for _ in 1:n_substeps
            rho_A = U_evol * rho_A * U_evol_adj
            rho_B = U_evol * rho_B * U_evol_adj
            rho_C = U_evol * rho_C * U_evol_adj
        end
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
    capacities, total_stm = calculate_stm_capacity(dict_B, inputs, MAX_DELAY, 50)
    plot_memory_capacity(capacities, MAX_DELAY, total_stm, Experiment_path)
end
run_task()