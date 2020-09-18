using Plots,FFTW,RecipesBase,Printf
export MvVmd,mvmd
"""
    MvVmd{T,S<:Int}

    MultiVariable Variational Mode Decompostion 
    Return by vmd function.

** Arguments **
- `signal::Array{T,2}`:signal
- `K::S`              :modes 
- `signal_d::Array{T,Val{K}}`:decomposed signals
- `freqs::Array{T,Val{K}}`:decomposed spectrums
- `samples::S`:signal sampled frequency
"""    
struct MvVmd{T,S<:Int}
    signal::Array{T,2}
    channels::S
    K::S
    signal_d::Array{T,3}
    mode::Array{Complex{T},3}
    omega::Array{T,1}
    sample_frequency::S
    MvVmd(signal,channels,K,decomp,mode,omega,sample_frequency)=new{eltype(signal),typeof(K)}(signal,channels,K,decomp,mode,omega,sample_frequency)
end


function Base.show(io::IO,v::MvVmd)
    println("Multiple Variable Variational Mode Decompostion")
    println("Channels : ",v.channels)
    println("Sample Frequency : $(v.sample_frequency) Hz")
    println("Base Frequency")
    for i in 1:v.K
        print("No $i component frequency = ")
        Printf.@printf("%5.3f Hz \n",v.omega[i])
    end
end
"""
 Multivariate Variational Mode Decomposition
    
     The function MVMD applies the "Multivariate Variational Mode Decomposition (MVMD)" algorithm to multivariate or multichannel data sets. 
     We have verified this code through simulations involving synthetic and real world data sets containing 2-16 channels. 
     However, there is no reason that it shouldn't work for data with more than 16 channels.
     
    
     Input and Parameters:
     ---------------------
     signal  - input multivariate signal that needs to be decomposed
     alpha   - the parameter that defines the bandwidth of extracted modes (low value of alpha yields higher bandwidth)
     tau     - time-step of the dual ascent ( pick 0 for noise-slack )
     K       - the number of modes to be recovered
     DC      - true if the first mode is put and kept at DC (0-freq)
     init    - 0 = all omegas start at 0
             - 1 = all omegas start uniformly distributed
             - 2 = all omegas initialized randomly
     tol     - tolerance value for convergence of ADMM
    
    
     Output:
     ----------------------
     u       - the collection of decomposed modes
     u_hat   - spectra of the modes
     omega   - estimated mode center-frequencies
    
    
     Syntax:
     -----------------------
 
     [u, u_hat, omega] = MVMD(X, alpha, tau, K, DC, init, tol)
       returns:
    			 a 3D matrix 'u(K,L,C)' containing K multivariate modes, each with 'C' number of channels and length 'L', that are 
                computed by applying the MVMD algorithm on the C-variate signal (time-series) X of length L.
        		 - To extract a particular channel 'c' corresponding to all extracted modes, you can use u_c = u(:,:,c).
    			 - To extract a particular mode 'k' corresponding to all channels, you can use u_k = u(k,:,:).
    			 - To extract a particular mode 'k' corresponding to the channel 'c', you can use u_kc = u(k,:,c).
    			 3D matrix 'u_hat(K,L,C)' containing K multivariate modes, each with 'C' number of channels and length 'L', that  
                are computed by applying the MVMD algorithm on the C-variate signal (time-series) X of length L.
    			 2D matrix 'omega(N,K)' estimated mode center frequencies
     Usage:
     -----------------------
     	Example 1: Mode Alignment on Synthetic Data
     	T = 1000; t = (1:T)/T;
     	f_channel1 = (cos(2*pi*2*t)) + (1/16*(cos(2*pi*36*t)));  Channel 1 contains 2 Hz and 36Hz tones
     	f_channel2 = (1/4*(cos(2*pi*24*t))) + (1/16*(cos(2*pi*36*t)));  Channel 2 contains 24 Hz and 36Hz tones
     	f = [f_channel1;f_channel2];  Making a bivariate signal
     	[u, u_hat, omega] = MVMD(f, 2000, 0, 3, 0, 1, 1e-7); 
     	Example 2: Real World Data (EEG Data)
     	load('EEG_data.mat');
     	[u, u_hat, omega] = MVMD(data, 2000, 0, 6, 0, 1, 1e-7);
     	Authors: Naveed ur Rehman and Hania Aftab
     	Contact Email: naveed.rehman@comsats.edu.pk
    
     	Acknowledgments: The MVMD code has been developed by modifying the univariate variational mode decomposition code that has 
                     been made public at the following link. We are also thankful to Dr. Maik Neukrich who helped us in developing a newer faster 
                     version of the code.  
                     https://www.mathworks.com/matlabcentral/fileexchange/44765-variational-mode-decomposition
                     by K. Dragomiretskiy, D. Zosso.
    		  
                     
    
     	Please cite the following papers if you use this code in your work:
       -----------------------------------------------------------------
     
      [1] N. Rehman, H. Aftab, Multivariate Variational Mode Decomposition, arXiv:1907.04509, 2019. 
      [2] K. Dragomiretskiy, D. Zosso, Variational Mode Decomposition, IEEE Transactions on Signal Processing, vol. 62, pp. 531-544, 2014. 
    ---------- Check for Input Signal              
     Check for getting number of channels from input signal
"""
function mvmd(signal::Array{Typ,2};alpha=2*size(signal,1), tau=0, K=3, DC=false, init=1, tol=1e-7, sample_frequency=100,iters=500) where Typ

    rows, cols = size(signal)
    @assert rows > cols "row should greater than column"
    
    #---------- Preparations
    #Sampling Frequency
    C = cols
    T = rows
    fs = 1/T
    T2 = T÷2
    #Mirroring
    f = [signal[T2:-1:1,:];signal;signal[T:-1:T2+1,:]]
    
    #Time Domain 0 to T (of mirrored signal)
    T = size(f,1)
    T2 = T ÷ 2
    t = collect(1:T)/T
    #frequencies
    freqs = t .- (0.5+1/T)
    #Construct and center f_hat
    f_hat_plus = fftshift(fft(f,1),1)
    f_hat_plus[1:T2,:] .= 0.0 + 0.0im
    
    #------------ Initialization
    #Maximum number of iterations 
    
    #For future generalizations: individual alpha for each mode
    Alpha = alpha*ones(K)
    #matrix keeping track of every iterant 
    u_hat_plus_m1 = zeros(Complex{Typ}, length(freqs), C, K)
    u_hat_plus_00 = zeros(Complex{Typ}, length(freqs), C, K)
    uDiff = zeros(Complex{Typ}, length(freqs), C, K)
    u_hat_plus = zeros(Complex{Typ},length(freqs), C, K)
    omega_plus = zeros(iters, K)
    #initialize omegas uniformly
    
    if init == 1
        omega_plus[1,:] .= (0.5/K)*(collect(1:K) .- 1.0)
    elseif inti == 2
        omega_plus[1,:] .= sort(exp(log(fs) + (log(0.5)-log(fs))*rand(K)))
    else
        omega_plus[1,:] .= 0
    end
    #if DC mode imposed, set its omega to 0
    if DC 
        omega_plus[1,1] = 0
    end
    #start with empty dual variables
    lambda_hat = zeros(Complex{Typ}, length(freqs), C, iters) 
    #other inits
    #update step
    n = 1  #loop counter
    sum_uk = zeros(Complex{Typ}, length(freqs), C)  #accumulator
    temp = zeros(Complex{Typ}, length(freqs), C)
    Δdiff = tol+eps()
    #--------------- Algorithm of MVMD
    while ( Δdiff > tol &&  n < iters )  #not converged and below iterations limit	
        #update modes
        for k ∈ 1:K
             #update mode accumulator
            if k > 1
                for i ∈ 1:length(freqs), c ∈ 1:C
                    sum_uk[i,c] += u_hat_plus[i,c,k-1] - u_hat_plus_00[i,c,k]
                end
            else
                for i ∈ 1:length(freqs), c ∈ 1:C
                    sum_uk[i,c] += u_hat_plus_00[i,c,K] - u_hat_plus_00[i,c,k]
                end
            end
            #update spectrum of mode through Wiener filter of residuals
            for i ∈ 1:length(freqs), c ∈ 1:C
                u_hat_plus[i,c,k] = (f_hat_plus[i,c] - sum_uk[i,c] - lambda_hat[i,c,n]/2)/(1+Alpha[k]*(freqs[i] - omega_plus[n,k])^2)
            end
            #update first omega if not held at 0
            if DC || (k > 1)
                #center frequencies
                numerator = reshape(freqs[T2+1:T],(1,length(freqs[T2+1:T])))*(abs.(u_hat_plus[T2+1:T,:, k]).^2)
                denominator = sum(abs.(u_hat_plus[T2+1:T,:,k]).^2)
                temp1 = sum(numerator)
                temp2 = sum(denominator)
                omega_plus[n+1,k] = temp1/temp2
            end
        end
        #Dual ascent
        temp .= reshape(sum(u_hat_plus, dims=3),(size(u_hat_plus,1),size(u_hat_plus,2)))
        for i ∈ 1:length(freqs), c ∈ 1:C
            lambda_hat[i,c,n+1] = lambda_hat[i,c,n] + tau*(temp[i,c] - f_hat_plus[i,c])
        end
        #loop counter
        n = n+1
        u_hat_plus_m1 .= u_hat_plus_00
        u_hat_plus_00 .= u_hat_plus
        #converged yet?
        uDiff .= u_hat_plus_00 - u_hat_plus_m1
        uDiff .= 1/T*(uDiff).*conj(uDiff)
        Δdiff = eps()+abs(sum(uDiff[:]))
    end
    #------ Post-processing and cleanup
    #discard empty space if converged early
    N = min(iters,n)
    omega = omega_plus[1:N,:]
    #Signal reconstruction
    u_hat = zeros(Complex{Typ},length(freqs), C, K)
    for c = 1:C
        u_hat[(T2+1):T,c,:] .= u_hat_plus[(T2+1):T,c,:]
        u_hat[(T2+1):-1:2,c,:] .= conj(u_hat_plus[(T2+1):T,c,:])
        u_hat[1,c,:] .= conj(u_hat[end,c,:])
    end
    u = zeros(Complex{Typ},length(freqs),C,K)
    for k = 1:K
        for c = 1:C
            u[:,c,k] .= real(ifft(ifftshift(u_hat[:,c,k])))
        end
    end
    #remove mirror part
    T4 = T ÷ 4
    u = u[T4+1:3*T4,:,:]
    u_hat = u_hat[T4+1:3*T4,:,:]
    #recompute spectrum
    #clear u_hat;

    for k = 1:K
        for c = 1:C
            u_hat[:,c,k] .= fftshift(fft(u[:,c,k]))
        end
    end
    println("iters = ",n)
    MvVmd(signal, C, K, u, u_hat, omega[end,:]*sample_frequency, sample_frequency)
end

"""
    plot(vmd::Vmd;k=1)
    
    visual the decomposed signals and spectrums

** Argument **
- `vmd::Vmd` : vmd
- `k::Int`   : 0-original signal 1-first component enforcement
"""
@recipe function f(v::MvVmd,::Val{:signal})
    
    T = size(v.signal,1)
    t = (1:T)./v.sample_frequency
    layout := (v.K+1, v.channels)
    for c ∈ 1:v.channels
        @series begin
            subplot := c
            xguide := "Time (s)"
            yguide := "c0"
            title  := "channel $c"
            t,v.signal[:,c]
        end
    end

    for k ∈ 1:v.K, c ∈ 1:v.channels
        @series begin
            subplot := k * v.channels + c
            xguide := "Time (s)"
            yguide := "c$k"
            label := "$(Printf.@sprintf("%5.3f",v.omega[k])) Hz"
            t,v.signal_d[:,c,k]
        end
    end
 end

 @recipe function f(v::MvVmd,::Val{:spectrum})
    
    T = size(v.signal,1)
    T2 = T ÷ 2
    t = collect(1:T)/T
    freqs = @. (t- 0.5 - 1 / T) * v.sample_frequency
    freqs = freqs[T2+1:end]
    layout := (v.K+1, v.channels)
    f = fftshift(fft(v.signal, 1), 1)
    for c ∈ 1:v.channels
        @series begin
            subplot := c
            xguide := "Time (s)"
            yguide := "c0 db"
            title  := "channel $c"
            freqs,20log10.(abs.(f[:,c][T2+1:T]))
        end
    end

    for k ∈ 1:v.K, c ∈ 1:v.channels
        @series begin
            subplot := k * v.channels + c
            xguide := "Time (s)"
            yguide := "c$k db"
            minorticks := false
            label := "$(Printf.@sprintf("%5.3f",v.omega[k])) Hz"
            freqs,20log10.(abs.(v.mode[:,c,k][T2+1:T]))
        end
    end
 end
 using Random
 T = 1000; t = (1:T)/T;
 f_channel1 = (cos.(2*pi*2*t)) + (1/16*(cos.(2*pi*60*t))) + 0.001 * randn(length(t));
 f_channel2 = (1/4*(cos.(2*pi*24*t))) + (1/16*(cos.(2*pi*36*t))) + 0.001 * randn(length(t));
 f = [f_channel1 f_channel2];
 v = mvmd(f;K = 4,DC = true, init = 1, sample_frequency = 1000)

 plot(v,Val(:signal))
 plot(v,Val(:spectrum))