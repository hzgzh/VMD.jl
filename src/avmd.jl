using Entropies
using StatsBase
using LinearAlgebra

function adaptive_vmd(signal; alpha=2000, tol=1e-7, max_iter=500,
                      K_max=10, threshold=1e-6,
                      m=3, τ=1)
    N = length(signal)
    scores = []
    K_values = []

    for K in 1:K_max
        v = vmd(signal; K=K, alpha=alpha)
        modes = v.u
        reconstructed_signal = sum(modes[:, k] for k in 1:K)
        E_loss = norm(signal - reconstructed_signal)^2 / norm(signal)^2

        if E_loss >= threshold
            continue
        end

        est = SymbolicPermutation(m=m, τ=τ)
        pe_values = [entropy(modes[:, k], est) for k in 1:K]
        PE_avg = mean(pe_values)

        kurt_values = [abs(StatsBase.kurtosis(modes[:, k])) for k in 1:K]
        avg_kurtosis = mean(kurt_values)

        score = PE_avg + avg_kurtosis + E_loss
        push!(scores, score)
        push!(K_values, K)
    end

    if isempty(scores)
        error("No suitable K found within given constraints")
    end
    min_index = argmin(scores)
    K_opt = K_values[min_index]

    v_opt = vmd(signal; K=K_opt, alpha=alpha)
    modes_opt = v_opt.u
    omega_opt = v_opt.w

    return modes_opt[:, 1:K_opt], omega_opt[1:K_opt], K_opt
end


function reconstruction_error(f::Vector{T}, u::Matrix{T}) where T
    reconstructed = sum(u[:, k] for k in 1:size(u)[2])
    err = norm(f - reconstructed)^2 / norm(f)^2
    return err
end

function band_widths(u::Matrix{T}, omega::Vector{T}, fs::Real=1.0) where T
    N = size(u)[1]
    freq = fftfreq(N) * fs
    band_widths_list = Vector{T}(undef, size(u)[2])
    for k in 1:size(u)[2]
        u_k_fft = fft(u[:, k])
        P_k_all = abs.(u_k_fft).^2
        mean_f_all = sum(freq .* P_k_all) / sum(P_k_all)
        var_f_all = sum((freq .- mean_f_all).^2 .* P_k_all) / sum(P_k_all)
        std_f_all = sqrt(var_f_all)
        band_widths_list[k] = std_f_all
    end
    return band_widths_list
end

function average_bandwidth(u::Matrix{T}, omega::Vector{T}, fs::Real=1.0) where T
    bw_vec = band_widths(u, omega, fs)
    return mean(bw_vec)
end

function average_correlation(u::Matrix{T}) where T
    num_modes = size(u)[2]
    if num_modes < 2
        return NaN
    end
    correlations_sum = zero(T)
    count_pairs = zero(Int)
    for i in 1:num_modes-1
        for j in i+1:num_modes
            cor_ij = abs(dot(u[:,i], u[:,j])) / (norm(u[:,i]) * norm(u[:,j]))
            correlations_sum += cor_ij
            count_pairs += 1
        end
    end
    return correlations_sum / count_pairs
end

function average_permutation_entropy(u::Matrix{T}; m=3, t=1) where T
    est = SymbolicPermutation(m=m, t=t)
    pe_values = [entropy(u[:,k], est) for k in 1:size(u)[2]]
    return mean(pe_values)
end

function mode_separation(omega::Vector{T}, bandwidths::Vector{T}; c=1) where T
    num_modes = length(omega)
    if num_modes < 2
        return NaN
    end
    separated_pairs_count = zero(Int)
    total_pairs_count = (num_modes*(num_modes-1)) ÷ 2
    for i in 1:num_modes-1
        for j in i+1:num_modes
            f_i = omega[i] / (2*pi)
            f_j = omega[j] / (2*pi)
            bw_i = bandwidths[i]
            bw_j = bandwidths[j]
            if abs(f_i - f_j) > c*(bw_i + bw_j)
                separated_pairs_count += 1
            end
        end
    end
    return separated_pairs_count / total_pairs_count
end

function cost_function(E_recon, B_avg, O_avg=nothing, PE_avg=nothing, MS=nothing;
                      lambda=1.0, MU=nothing, NU=nothing, GAMMA=nothing)
    cost = E_recon + lambda * B_avg
    if !isnothing(O_avg) && !isnothing(MU)
        cost += MU * O_avg
    end
    if !isnothing(PE_avg) && !isnothing(NU)
        cost += NU * PE_avg
    end
    if !isnothing(MS) && !isnothing(GAMMA)
        cost -= GAMMA * MS
    end
    return cost
end

function optimize_vmd_params(signal; K_range=2:8, alpha_range=500:500:5000, fs=1.0,
                            lambda=1.0, MU=nothing, NU=nothing, GAMMA=nothing,
                            m=3, t=1, c=1)
    best_cost = Inf
    best_K = nothing
    best_alpha = nothing
    best_modes = nothing
    best_omega = nothing
    
    for K in K_range
        for alpha in alpha_range
            try
                vmd_result = vmd(signal; K=K, alpha=alpha)
                u = vmd_result.u
                omega = vmd_result.w
                
                E_recon = reconstruction_error(signal, u)
                bw_vec = band_widths(u, omega, fs)
                B_avg = mean(bw_vec)
                O_avg = average_correlation(u)
                PE_avg = average_permutation_entropy(u, m=m, t=t)
                MS = mode_separation(omega, bw_vec, c=c)
                
                cost = cost_function(E_recon, B_avg, O_avg=O_avg, PE_avg=PE_avg, MS=MS,
                                    lambda=lambda, MU=MU, NU=NU, GAMMA=GAMMA)
                
                if cost < best_cost
                    best_cost = cost
                    best_K = K
                    best_alpha = alpha
                    best_modes = u
                    best_omega = omega
                end
                
            catch e
                @error "Error with K=$K and alpha=$alpha: $e"
            end
        end
    end
    
    return best_K, best_alpha, best_modes, best_omega, best_cost
end
