using VMD
using Test
using Random


@testset "VMD.jl" begin
    # Write your tests here.
    T = 1000
    t = (1:T) / T
    sample_frequency = 1000

# center frequencies of components
    f_1 = 2
    f_2 = 24
    f_3 = 288
    

# modes
    v_1 = @. cos(2 * pi * f_1 * t)
    v_2 = @. 1 / 4 * (cos(2 * pi * f_2 * t))
    v_3 = @. 1 / 16 * (cos(2 * pi * f_3 * t))

# composite signal, including noise
    f = v_1 + v_2 + v_3 + 0.1 * randn(length(v_1))

# some sample parameters for VMD
    alpha = 2000       # moderate bandwidth constraint
    tau = 0            # noise-tolerance (no strict fidelity enforcement)
    K = 4             # 3 modes + noise
    DC = false             # no DC part imposed
    init = 1           # initialize omegas uniformly
    tol = 1e-7


    v = vmd(f ; alpha=alpha,tau=tau,K=K,DC=false,init=1,tol=tol,sample_frequency=sample_frequency)

    @test isapprox(v.omega[1] , f_1,rtol=0.01)
    @test isapprox(v.omega[2] , f_2,rtol=0.01)
    @test isapprox(v.omega[3] , f_3,rtol=0.01)
end
