using LinearAlgebra

# --- IMPORTACIÓN DE MÓDULOS ---
include("src/operator_terms/pauli.jl")      
include("src/operator_terms/hamiltonian.jl")
include("src/utils/dynamics.jl")
include("src/utils/injection.jl")           
include("src/visualization/plot_spin_dynamics.jl") 

function run_online_protocol()
    # 1. Configuración
    n_qubits = 6
    J = 1.0     
    h = 0.5     # Campo fuerte para mezclar
    dt = 0.05   
    
    # Duración de la evolución entre cada input
    steps_per_input = 10000 
    
    println("Configurando QRC Online (Schrödinger Forward)...")
    
    # 2. Hamiltoniano (Truco: Signo negativo para evolucionar estados rho)
    H_ising = build_ising_hamiltonian(n_qubits, J, h)
    H_schrodinger = Operator()
    for (p, c) in H_ising
        H_schrodinger[p] = -c 
    end
    
    # 3. Inputs: Onda cuadrada (+0.8, -0.8...)
    n_inputs = 300
    inputs = [isodd(i) ? 0.9 : -0.9 for i in 1:n_inputs]
    
    # 4. Estado Inicial (rho)
    # Empezamos con el sistema "vacío" (Identidad global ~ Temperatura infinita)
    # y le inyectamos el primer input al Qubit 0.
    rho = Operator()
    # El término identidad global siempre está implícito en Pauli, 
    # así que empezamos definiendo explícitamente el estado del Qubit 0.
    rho[PauliString(0, 1)] = inputs[1] + 0.0im 
    
    # Matrices para guardar datos
    times = Float64[]
    # Guardaremos Z0, Z1... y Z5
    results = zeros(Float64, n_inputs, n_qubits)
    
    println("Procesando $n_inputs inputs...")
    
    for k in 1:n_inputs
        input_u = inputs[k]
        
        # A) INYECCIÓN (ERASE & WRITE)
        # Excepto en el paso 1 (que ya inicializamos arriba), reseteamos el Q0
        if k > 1
            rho = inject_input_schrodinger(rho, 0, input_u)
        end
        
        # B) EVOLUCIÓN (Reservoir Mixing)
        for s in 1:steps_per_input
            rho = step_rk4(rho, H_schrodinger, dt)
        end
        
        # C) MEDICIÓN (Features)
        push!(times, k) # Eje X = Número de input
        
        for q in 0:(n_qubits-1)
            # Medir Z en qubit q
            val = get(rho, PauliString(0, 1 << q), 0.0im)
            results[k, q+1] = real(val)
        end
    end
    
    println("Protocolo finalizado.")
    
    # Visualizar (Usamos n_singles=n_qubits, n_corrs=0 para ver solo los Z individuales)
    labels = ["Z$i" for i in 0:(n_qubits-1)]
    plot_dynamics_custom(times, results, labels, n_qubits, 0)
    
    # Imprimir últimos valores para verificar
    println("Valores finales (Input fue $(inputs[end])): ")
    println(round.(results[end, :], digits=4))
end

run_online_protocol()