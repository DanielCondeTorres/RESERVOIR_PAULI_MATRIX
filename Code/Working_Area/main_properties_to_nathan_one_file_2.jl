using LinearAlgebra
using Plots
using Random
using JLD2
using FileIO
using Statistics

# ==============================================================================
# 1. MOTOR MATRICIAL Y OPERADORES
# ==============================================================================
const SCRIPT_DIR = @__DIR__
const I2 = ComplexF64[1 0; 0 1]
const X  = ComplexF64[0 1; 1 0]
const Y  = ComplexF64[0 -im; im 0]
const Z  = ComplexF64[1 0; 0 -1]
const Zero_Proj = ComplexF64[1 0; 0 0] # |0><0|
const Flip_Op   = ComplexF64[0 1; 0 0] # |0><1|

# Construcción de operadores (Little-Endian: Q0 a la derecha)
function build_op(n_qubits, target_idx, gate)
    res = (target_idx == 0) ? gate : I2
    for i in 1:(n_qubits - 1)
        res = kron((i == target_idx) ? gate : I2, res)
    end
    return res
end

function build_nathan_hamiltonian(n_qubits::Int, h::Float64; seed::Int=12345)
    Random.seed!(seed)
    dim = 2^n_qubits
    H = zeros(ComplexF64, dim, dim)
    for j in 0:(n_qubits-1)
        H -= h * build_op(n_qubits, j, Z)
    end
    for i in 0:(n_qubits-1), j in (i+1):(n_qubits-1)
        J_ij = (rand() - 0.5)
        H += J_ij * (build_op(n_qubits, i, X) * build_op(n_qubits, j, X))
    end
    return H
end

function inject_state_EraseWrite(rho::Matrix{ComplexF64}, qubit_idx::Int, 
    rz::Float64, rx::Float64=0.0, ry::Float64=0.0)
    N = Int(log2(size(rho, 1)))
    M0 = build_op(N, qubit_idx, Zero_Proj)
    M1 = build_op(N, qubit_idx, Flip_Op)
    rho_reset = M0 * rho * M0' + M1 * rho * M1'
    norm_v = sqrt(rx^2 + ry^2 + rz^2)
    if norm_v < 1e-9; return rho_reset; end
    theta, phi = acos(clamp(rz/norm_v, -1.0, 1.0)), atan(ry, rx)
    U_rot = [exp(-im*phi/2) 0; 0 exp(im*phi/2)] * [cos(theta/2) -sin(theta/2); sin(theta/2) cos(theta/2)]
    return build_op(N, qubit_idx, U_rot) * rho_reset * build_op(N, qubit_idx, U_rot)'
end

# ==============================================================================
# 2. TUS FUNCIONES DE CAPACIDAD (INTEGRADAS)
# ==============================================================================

function train_reservoir(observables::Matrix{Float64}, targets::Vector{Float64}, washout::Int=20)
    X = observables[washout+1:end, :]
    y = targets[washout+1:end]
    X_bias = hcat(X, ones(size(X, 1)))
    weights = X_bias \ y
    return weights, X_bias * weights
end

function calculate_capacity(target::Vector{Float64}, prediction::Vector{Float64})
    n = min(length(target), length(prediction))
    y, y_pred = target[end-n+1:end], prediction[end-n+1:end]
    cv, v_y, v_pred = cov(y, y_pred), var(y), var(y_pred)
    return (v_y < 1e-12 || v_pred < 1e-12) ? 0.0 : (cv^2) / (v_y * v_pred)
end

function calculate_stm_capacity(reservoir_data::Dict{String, Vector{Float64}}, original_inputs::Vector, max_delay::Int=15, washout::Int=50)
    targets_float = Float64.(original_inputs)
    steps = length(targets_float)
    feature_keys = sort(collect(keys(reservoir_data)))
    X_reservoir = hcat([reservoir_data[k] for k in feature_keys]...)
    
    capacities = Float64[]
    for tau in 0:max_delay
        if tau == 0
            y_target, X_feats = targets_float, X_reservoir
        else
            y_target, X_feats = targets_float[1:end-tau], X_reservoir[1+tau:end, :]
        end
        if length(y_target) <= washout; push!(capacities, 0.0); continue; end
        
        X_clean, y_clean = X_feats[washout+1:end, :], y_target[washout+1:end]
        n_split = floor(Int, length(y_clean) * 0.8)
        if n_split < 10; push!(capacities, 0.0); continue; end
        
        weights, _ = train_reservoir(X_clean[1:n_split, :], y_clean[1:n_split], 0)
        y_pred_test = hcat(X_clean[n_split+1:end, :], ones(size(X_clean, 1)-n_split)) * weights
        push!(capacities, calculate_capacity(y_clean[n_split+1:end], y_pred_test))
    end
    return capacities, sum(capacities)
end

function plot_memory_capacity(capacities::Vector{Float64}, max_delay::Int, total_stm::Float64, save_dir::String)
    taus = 0:length(capacities)-1
    p = bar(taus, capacities, label="Memory Capacity", xlabel="Delay (τ)", ylabel="Capacity C(τ)",
        title="STM Capacity (Total = $(round(total_stm, digits=2)))", color=:skyblue, ylims=(0, 1.1), alpha=0.8)
    plot!(p, taus, capacities, color=:blue, label="", lw=2)
    savefig(p, joinpath(save_dir, "STM_Capacity_Plot.png"))
end

# ==============================================================================
# 3. REPORTES Y SIMULACIÓN
# ==============================================================================

function plot_and_save_validation_full(dict_A, dict_B, N, steps, save_dir; label_A="Tray. A", label_B="Tray. B")
    basis = any(k->contains(k,"X"), keys(dict_A)) ? "X" : (any(k->contains(k,"Y"), keys(dict_A)) ? "Y" : "Z")
    mkpath(joinpath(save_dir, "Individual_$basis")); mkpath(joinpath(save_dir, "Pairs_$basis$basis")); mkpath(joinpath(save_dir, "Triads_$basis$basis$basis"))
    for lbl in keys(dict_A)
        n_ops = count(c -> string(c) == basis, lbl)
        idxs = [i for (i, c) in enumerate(lbl) if string(c) == basis]
        tag = basis * join(idxs)
        p = plot(dict_A[lbl], label=label_A, c=:blue, title=tag, lw=1.2)
        plot!(p, dict_B[lbl], label=label_B, c=:red, ls=:dash, lw=1.2)
        ylims!(p, -1.1, 1.1)
        subf = (n_ops == 1) ? "Individual_$basis" : (n_ops == 2 ? "Pairs_$basis$basis" : "Triads_$basis$basis$basis")
        savefig(p, joinpath(save_dir, subf, "Plot_$tag.png"))
    end
end

function run_task()
    N, steps, T_evol, h_val = 6, 200, 100.0, 1.0 # Aumenté pasos para que el split 80/20 tenga sentido
    Experiment_path = joinpath(SCRIPT_DIR, "Experiment_Nathan_Final")
    mkpath(Experiment_path)

    H = build_nathan_hamiltonian(N, h_val)
    U_evol = exp(-im * H * T_evol); U_adj = adjoint(U_evol)

    rho_A, rho_B = zeros(ComplexF64, 2^N, 2^N), zeros(ComplexF64, 2^N, 2^N)
    rho_A[1,1] = 1.0; rho_B[end,end] = 1.0
    rho_C = copy(rho_A)

    Random.seed!(1234); inputs = rand(steps)
    inputs_C = copy(inputs); if steps > 30; inputs_C[31:end] = rand(steps-30); end

    obs_matrices = Dict{String, Matrix{ComplexF64}}()
    for (b_char, b_op) in zip(['X','Y','Z'], [X, Y, Z])
        for i in 0:N-1; obs_matrices[join([m==i ? b_char : '1' for m in 0:N-1])] = build_op(N, i, b_op); end
        for i in 0:N-1, j in i+1:N-1; obs_matrices[join([m==i||m==j ? b_char : '1' for m in 0:N-1])] = build_op(N, i, b_op) * build_op(N, j, b_op); end
        for i in 0:N-1, j in i+1:N-1, k in j+1:N-1; obs_matrices[join([m==i||m==j||m==k ? b_char : '1' for m in 0:N-1])] = build_op(N, i, b_op) * build_op(N, j, b_op) * build_op(N, k, b_op); end
    end

    dict_A = Dict(l => zeros(steps) for l in keys(obs_matrices))
    dict_B = Dict(l => zeros(steps) for l in keys(obs_matrices))
    dict_C = Dict(l => zeros(steps) for l in keys(obs_matrices))

    println("🏃 Simulando...")
    for s in 1:steps
        rz, rz_C = 1.0 - 2.0*inputs[s], 1.0 - 2.0*inputs_C[s]
        rho_A = inject_state_EraseWrite(rho_A, 0, rz)
        rho_B = inject_state_EraseWrite(rho_B, 0, rz)
        rho_C = inject_state_EraseWrite(rho_C, 0, rz_C)
        rho_A, rho_B, rho_C = U_evol*rho_A*U_adj, U_evol*rho_B*U_adj, U_evol*rho_C*U_adj
        for (lbl, op) in obs_matrices
            dict_A[lbl][s], dict_B[lbl][s], dict_C[lbl][s] = real(tr(rho_A*op)), real(tr(rho_B*op)), real(tr(rho_C*op))
        end
    end

    for b in ["X", "Y", "Z"]
        b_dir = joinpath(Experiment_path, "Basis_$b"); mkpath(b_dir)
        plot_and_save_validation_full(filter(p->contains(p.first, b), dict_A), filter(p->contains(p.first, b), dict_B), N, steps, joinpath(b_dir, "Echo"))
        plot_and_save_validation_full(filter(p->contains(p.first, b), dict_A), filter(p->contains(p.first, b), dict_C), N, steps, joinpath(b_dir, "Separability"), label_B="Tray. C")
    end

    caps, tot = calculate_stm_capacity(dict_A, inputs, 10, 50)
    plot_memory_capacity(caps, 10, tot, Experiment_path)
    println("✅ Finalizado. STM Total: $(round(tot, digits=3))")
end

run_task()