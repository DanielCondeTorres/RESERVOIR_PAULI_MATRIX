using LinearAlgebra

# --- IMPORTACIÓN DE MÓDULOS ---
include("src/operator_terms/pauli.jl")      
include("src/operator_terms/hamiltonian.jl")
include("src/utils/dynamics.jl")
include("src/utils/injection.jl")           
include("src/visualization/plot_spin_dynamics.jl") 
function run_online_protocol()
    # --- CONFIGURACIÓN ---
    n_qubits = 6
    J = 1.0     
    h = 0.5     # Campo transversal
    dt = 0.05   
    
    # Parámetros de Ruido (CRÍTICO)
    gamma_dephasing = 0.005 # Tasa de decaimiento por paso de tiempo
    
    steps_per_input = 200 # Reducido de 10000 (demasiado largo) a 200 para ver dinámica
    n_inputs = 100
    
    println("Configurando QRC Online con inject_state unificado...")
    
    # --- HAMILTONIANO ---
    # Usamos signo negativo porque dRho/dt = -i[H, Rho]. 
    # En RK4 simple a veces se integra H directo, aquí asumimos tu convención.
    H_ising = build_ising_hamiltonian(n_qubits, J, h)
    H_schrodinger = Operator()
    for (p, c) in H_ising
        H_schrodinger[p] = -c 
    end
    
    # --- INPUTS ---
    # Onda cuadrada (+0.8, -0.8...)
    inputs = [isodd(i) ? 0.9 : -0.9 for i in 1:n_inputs]
    
    # --- ESTADO INICIAL ---
    rho = Operator()
    # Inicialización limpia usando la NUEVA FUNCIÓN
    # Inyectamos el primer input en el Qubit 0 (índice 0)
    rho = inject_state(rho, 0, inputs[1]) 
    
    # Matrices para guardar datos
    times = Float64[]
    results = zeros(Float64, n_inputs, n_qubits)
    
    println("Procesando $n_inputs inputs...")
    
    for k in 1:n_inputs
        input_u = inputs[k]
        
        # A) INYECCIÓN (ERASE & WRITE)
        # Usamos inject_state. Automáticamente asigna input_u a rz.
        if k > 1
            rho = inject_state(rho, 0, input_u)
        end
        
        # B) EVOLUCIÓN + RUIDO
        for s in 1:steps_per_input
            # 1. Evolución Unitaria (Hamiltoniano)
            rho = step_rk4(rho, H_schrodinger, dt)
            
            # 2. Evolución No-Unitaria (Dephasing)
            # Esto es lo que permite olvidar el input anterior (Z decay)
            rho = apply_global_dephasing(rho, gamma_dephasing)
        end
        
        # C) MEDICIÓN
        push!(times, k) 
        
        for q in 0:(n_qubits-1)
            # Medir magnetización Z en qubit q
            # Nota: 1 << q funciona para índices 0, 1, 2...
            val = get(rho, PauliString(0, 1 << q), 0.0im)
            results[k, q+1] = real(val)
        end
    end
    
    println("Protocolo finalizado.")
    
    # Imprimir últimos valores
    println("Valores finales Z (Input actual: $(inputs[end])): ")
    println(round.(results[end, :], digits=4))
    
    # Descomentar si tienes la función de plot cargada
    labels = ["Z$i" for i in 0:(n_qubits-1)]
    plot_expectation_evolution(times, results, labels, n_qubits, 0)
end

run_online_protocol()