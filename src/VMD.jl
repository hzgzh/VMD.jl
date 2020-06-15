module VMD

using FFTW
using Plots
using Printf

export Vmd,vmd,plot,compare,n_component,n_mode

"""
    Vmd{T,S<:Int}

    Variational Mode Decompostion 
    Return by vmd function.

** Arguments **
- `signal::Array{T,1}`:signal
- `K::S`              :modes 
- `signal_d::Array{T,Val{K}}`:decomposed signals
- `freqs::Array{T,Val{K}}`:decomposed spectrums
- `samples::S`:signal sampled frequency

"""    
struct Vmd{T,S<:Int}
    signal::Array{T,1}
    K::S
    signal_d::Array{T,2}
    mode::Array{Complex{T},2}
    omega::Array{T,1}
    sample_frequency::S
    Vmd(signal,K,decomp,mode,omega,sample_frequency)=new{eltype(signal),typeof(K)}(signal,K,decomp,mode,omega,sample_frequency)
end

function Base.show(io::IO,v::Vmd)
    println("Variational Mode Decompostion")
    println("Sample Frequency : $(v.sample_frequency) Hz")
    println("Base Frequency")
    for i in 1:v.K
        print("No $i component frequency = ")
        Printf.@printf("%5.3f Hz \n",v.omega[i])
    end
end

"""
    n_mode(v::Vmd,k)
    return the k mode frequency
"""
n_mode(v::Vmd,k) = v.omega[k]

"""
    n_component(v::Vmd,k)

    return No k component decomposed signal
"""
n_component(v::Vmd,k) = v.signal_d[:,k] 

"""
    vmd(signal;alpha=100, tau=0, K=3, DC=false, init=1, tol=1e-6,samples=100)

    Variational Mode Decomposition
    Return Vmd

** Argument **

 - `signal`: the time domain signal (1D) to be decomposed
 - `alpha` : the balancing parameter of the data-fidelity constraint
 - `tau`   : time-step of the dual ascent ( pick 0 for noise-slack )
 - `K`     : the number of modes to be recovered
 - `DC`    : true if the first mode is put and kept at DC (0-freq)
 - `init`  : 0 = all omegas start at 0
             1 = all omegas start uniformly distributed
             2 = all omegas initialized randomly
 - `tol`   : tolerance of convergence criterion; typically around 1e-6
 - `sample_frequency` : samples frequency(eg:100/s)

 ** Example **
 ```julia
T = 1000;
t = (1:T)/T;
sample_frequency = 1000;

# center frequencies of components
f_1 = 2;
f_2 = 24;
f_3 = 288;

# modes
v_1 = @. cos(2*pi*f_1*t);
v_2 = @. 1/4*(cos(2*pi*f_2*t));
v_3 = @. 1/16*(cos(2*pi*f_3*t));

# composite signal, including noise
f = v_1 + v_2 + v_3 + 0.1*randn(length(v_1));

# some sample parameters for VMD
alpha = 2000;       # moderate bandwidth constraint
tau = 0;            # noise-tolerance (no strict fidelity enforcement)
K = 3;              # 3 modes
DC = false;             # no DC part imposed
init = 1;           # initialize omegas uniformly
tol = 1e-7;


v = vmd(f ; alpha = alpha,tau = tau,K = K,DC = false,init = 1,tol = tol,sample_frequency = sample_frequency)

# plot original signal and spectrum
plot(v;k = 0)

# plot first decomposed component and spectrum
plot(v;k = 1)
```
"""
function vmd(signal::Array{Typ,1};alpha=2*length(signal), tau=0, K=3, DC=false, init=1, tol=1e-6,sample_frequency=100) where Typ
    # ---------- Preparations

    # Period and sampling frequency of input signal
   
   
    T =  length(signal)
    fs = 1/T
    T2 = Int(T / 2)
     # extend the signal by mirroring
    f = [signal[T2:-1:1];signal;signal[T:-1:T2 + 1]]

    # Time Domain 0 to T (of mirrored signal)
    T = length(f)
    T2 = Int(T / 2)
    t = collect(1:T) / T

    # Spectral Domain discretization
    freqs = @. t - 0.5 - 1 / T

    # Maximum number of iterations (if not converged yet, then it won't anyway)
    N = 500

    # For future generalizations: individual alpha for each mode
    Alpha = alpha * ones(K)

    # Construct and center f_hat
    f_hat = fftshift((fft(f)))
    f_hat_plus = f_hat
    fill!(f_hat_plus[1:T2] ,0.0 + 0.0im)

    # matrix keeping track of every iterant // could be discarded for mem
    u_hat_plus = zeros(Complex{Typ},T, K,2)

    # Initialization of omega_k
    omega_plus = zeros(K,2)
    if init == 1
        for i = 1:K
            omega_plus[i,1] = (0.5 / K) * (i - 1)
        end
    elseif init == 2
        omega_plus[:,1] = sort(exp.(log(fs) .+ (log(0.5) - log(fs)) .* rand(K)))
    else
        omega_plus[:,1] .= 0
    end

    # if DC mode imposed, set its omega to 0
    if DC
        omega_plus[1, 1] = 0
    end

    # start with empty dual variables
    lambda_hat = zeros(Complex{Typ}, T,2)
    tmp = zeros(Complex{Float64}, T)
    # other inits
    uDiff = tol + eps() # update step
    n = 1 # loop counter
    sum_uk = zeros(Complex{Typ}, T) # accumulator

    # ----------- Main loop for iterative updates

    while abs(uDiff) > tol && n < N  # not converged and below iterations limit
        
        # update first mode accumulator
        k = 1
        
        @inbounds for i in 1:T
            sum_uk[i] +=  u_hat_plus[i,K,1] - u_hat_plus[i, 1 ,1]
                # update spectrum of first mode through Wiener filter of residuals
        end
        @inbounds for i in 1:T
            u_hat_plus[i,k,2]=(f_hat_plus[i]-sum_uk[i]-lambda_hat[i,1]/2.0)/(1.0 + Alpha[k]*(freqs[i]-omega_plus[k,1])^2)
        end
        
        a=zeros(Complex{Typ},2)
        # update first omega if not held at 0
        if ~DC
            @inbounds for i in T2+1:T
            a[1]+=freqs[i]*abs(u_hat_plus[i, k,2])^2
            a[2]+=abs(u_hat_plus[i, k,2])^2
            end
            omega_plus[k,2] =  a[1] / a[2]
        end

        # update of any other mode
        for k = 2:K

            # accumulator
            @inbounds for i in 1:T
                sum_uk[i] += u_hat_plus[i, k - 1,2] - u_hat_plus[i, k , 1]

            # mode spectrum
            end

            @inbounds for i in 1:T
                u_hat_plus[i, k , 2] = (f_hat_plus[i] - sum_uk[i] - lambda_hat[i , 1] / 2) / (1 +
                                     Alpha[k] * (freqs[i] - omega_plus[k , 1])^2)
            end
            # center frequencies
            a=zeros(Complex{Typ},2)
            @inbounds for i in T2+1:T
                a[1]+=freqs[i]*abs(u_hat_plus[i, k,2])^2
                a[2]+=abs(u_hat_plus[i, k,2])^2
            end
            omega_plus[k,2] =  a[1] / a[2]
        end

        # Dual ascentu
        @inbounds for i in 1:T,j=1:K
            tmp[i]+=u_hat_plus[i,j,2]
        end

        @inbounds for i in 1:T
             lambda_hat[i,2] = lambda_hat[i,1] + tau * (tmp[i] - f_hat_plus[i])
        end

        fill!(tmp,0.0+0.0im)
        # loop counter
        n = n + 1
        

        # converged yet?
        uDiff = eps()
        
        @inbounds for i in 1:K , j in 1:T
            uDiff += 1 / T * real((u_hat_plus[j, i ,2] - u_hat_plus[j, i ,1])*conj(u_hat_plus[j, i ,2] - u_hat_plus[j, i ,1]))
        end
        

        # prepare for iter
        @inbounds for i in 1:T
            for j in 1:K
                u_hat_plus[i,j,1] = u_hat_plus[i,j,2]
            end
            lambda_hat[i,1] = lambda_hat[i,2]
        end
        for i in 1:K
            omega_plus[i,1] = omega_plus[i,2]
        end
    end

    
    # ------ Postprocessing and cleanup
    

    # discard empty space if converged early
    
    omega = omega_plus[:,2]

    # Signal reconstruction
    u_hat = zeros(Complex{Typ}, T, K)
    
    u_hat[T2 + 1:T, :] = u_hat_plus[T2 + 1:T, :,2]
    u_hat[T2 + 1:-1:2, :] = conj(u_hat_plus[T2 + 1:T, :,2])
    u_hat[1, :] = conj(u_hat[end, :])

    u = zeros(T,K)
    
    for k = 1:K
        u[:,k] = real(ifft(ifftshift(u_hat[:, k])))
    end
    T4 = convert(typeof(T),T/4)
    # remove mirror part
    u = u[T4 + 1:3 * T4,:]

    # recompute spectrum
    # clear u_hat;
    u_hat = zeros(Complex{Typ}, size(u)[1], K)
    for k = 1:K
        u_hat[:, k] = fftshift(fft(u[:,k]))
    end
    println("--iteration times $n -- error $(abs(uDiff))")
    idx = sortperm(omega)
    omega = omega[idx]*sample_frequency
    u = u[:,idx]
    u_hat = u_hat[:,idx]

    Vmd(signal,K,u,u_hat,omega,sample_frequency)
end

"""
    plot(vmd::Vmd;k=1)
    
    visual the decomposed signals and spectrums

** Argument **
- `vmd::Vmd` : vmd
- `k::Int`   : 0-original signal 1-first component enforcement
"""
function Plots.plot(v::Vmd;k=1)
    @assert k<=v.K error("can't great than $(vmd.K)")

    T = length(v.signal)
    T2 = Int(T/2)
    t = collect(1:T)/T
    freqs = @. (t - 0.5 - 1 / T)*v.sample_frequency
    t = t*T/v.sample_frequency
    if k == 0
        f = fftshift(fft(v.signal))
        p1 = Plots.plot(t,v.signal,title="origin signal ",xlabel="Time (s)",ylabel = "y")
        p2 = Plots.plot(freqs[1:T2],20log.(abs.(f[T2+1:end])),title = "Original spectral",xlabel = "Freq Hz",ylabel = "db")
    else
        p1=Plots.plot(t,v.signal_d[:,k],title="Reconstructed mode $k",xlabel="Time (s)",ylabel = "y")
        p2=Plots.plot(freqs[T2+1:T],20log.(abs.(v.mode[T2+1:end,k])),title = "Spectral decomposition",
            xlabel = "Freq Hz",ylabel = "db",label = "$(Printf.@sprintf("%5.3f",v.omega[k])) Hz")
    end
    Plots.plot(p1,p2,layout = (2,1),size = (640,480))
end

"""
    compare(v::Vmd)
    
    visual compare origin signal and sum of decomposed signal 
"""
function compare(v::Vmd)
    T = length(v.signal)
    t = collect(1:T)./v.sample_frequency
    Plots.plot(t,v.signal,label = "origin signal")
    Plots.plot!(t,[sum(v.signal_d[i,:]) for i in 1:length(v.signal)],label = "reconstructed signal",xlabel="Time s")
end

end
