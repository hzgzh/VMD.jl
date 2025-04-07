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


function optimize_vmd_params(signal; K_range=2:8, alpha_range=500:500:5000, fs=1.0, lambda=1.0)
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
                
                # 计算重构误差
                E_recon = norm(signal - sum(u, dims=2))^2 / norm(signal)^2
                
                # 计算平均带宽（简化版）
                bw_vec = [sqrt(sum((fftfreq(length(signal))*fs .- omega[k]).^2 .* abs.(fft(u[:,k])).^2) / sum(abs.(fft(u[:,k])).^2)) for k in 1:K]
                B_avg = mean(bw_vec)
                
                # 成本函数
                cost = E_recon + lambda * B_avg
                
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
