using LinearAlgebra
using Statistics
using Plots
"""
Calculate the short-term memory (STM) capacity of a reservoir.
# Arguments
- `states_matrix`: Matrix of size (N_steps, N_observables) containing the reservoir states over time.
- `inputs`: Vector of size (N_steps) containing the input signal injected into the reservoir.
- `tau_max`: Maximum delay to calculate the capacity for.
# Returns
- `capacities`: Vector of size `tau_max` containing the memory capacity C(tau).
"""
function calculate_stm_capacity(states_matrix::Matrix{Float64}, inputs::Vector{Float64}, tau_max::Int=10; washout::Int=0)
    N_steps = length(inputs)
    capacities = zeros(Float64, tau_max)
    
    # Use data after the washout period to avoid transient effects
    valid_steps = N_steps - washout
    if valid_steps <= tau_max
        error("Not enough steps after washout to calculate capacity for tau_max = \$tau_max")
    end
    for tau in 1:tau_max
        # y is the target (delayed input)
        # We want to predict y_k = inputs[k - tau] using states_matrix[k, :]
        # We consider indices k from (tau + washout + 1) to N_steps
        
        start_idx = max(tau + 1, washout + 1)
        
        Y_target = inputs[start_idx-tau : N_steps-tau]
        X_train = states_matrix[start_idx : N_steps, :]
        
        # Add a bias column to X_train
        X_train_bias = hcat(ones(size(X_train, 1)), X_train)
        
        # Linear regression: W = (X^T X)^(-1) X^T Y
        # Compute pseudo-inverse to handle rank-deficient matrices
        W = pinv(X_train_bias) * Y_target
        
        # Predict
        Y_pred = X_train_bias * W
        
        # Calculate C(tau) using squared correlation coefficient (R^2 equivalent if zero mean, or squared Pearson correlation)
        cov_val = cov(Y_target, Y_pred)
        var_target = var(Y_target)
        var_pred = var(Y_pred)
        
        if var_target > 1e-12 && var_pred > 1e-12
            capacities[tau] = (cov_val^2) / (var_target * var_pred)
        else
            capacities[tau] = 0.0
        end
    end
    
    return capacities
end
function plot_stm_capacity(capacities::Vector{Float64}, out_path::String)
    tau_vals = 1:length(capacities)
    p = plot(tau_vals, capacities, 
             marker=:circle, 
             label="Capacity", 
             xlabel="\\tau", 
             ylabel="Capacity", 
             title="Short-term Memory Capacity",
             lw=2, color=:dodgerblue,
             legend=false, grid=true)
    
    savefig(p, out_path)
    println("✅ Gráfica de Capacidad STM guardada en: \$out_path")
    return p
end
