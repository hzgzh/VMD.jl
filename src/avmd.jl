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
