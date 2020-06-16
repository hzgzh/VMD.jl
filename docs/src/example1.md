## Example1

The first signal is a composition of three
simple components, namely a general linear trend and two dif-
ferent harmonics:

$$f_{\mathrm{Sig} 1}(t)=6 t+\cos (8 \pi t)+0.5 \cos (40 \pi t)$$

The signal, its three constituent modes, and the composite
Fourier spectrum are shown in Fig. 9. The main challenge
of this signal is the linear growth term, whose higher order
harmonics spread over the whole spectrum.

The recovered VMD modes constitute a nice partition of the
input spectrum, with each mode being clearly dominant around
its respective center frequency. The three modes in time domain
show nice separation into three distinct signals of characteristic
oscillations. Compared to the EMD results, the proposed VMD
method performs clearly better, with less spurious oscillations
in the trend and mid-frequency signal.

```@repl
using VMD,Random,Plots
T = 1000;
t = (1:T)/T;
sample_frequency = 1000;

# modes
v_1 = @. 6t;
v_2 = @. cos(8π*t)
v_3 = @. 0.5cos(40π*t) 

# composite signal, including noise
f = v_1 + v_2 + v_3 + 0.1*randn(length(v_1));

# some sample parameters for VMD
alpha = 2000;       # moderate bandwidth constraint
tau = 0;            # noise-tolerance (no strict fidelity enforcement)
K = 3;              # 3 modes
DC = false;             # no DC part imposed
init = 0;           # initialize omegas uniformly
tol = 1e-7;


v = vmd(f ; alpha = alpha,tau = tau,K = K,DC = false,init = init,tol = tol,sample_frequency = sample_frequency)


p0 = VMD.plot(v,k=0)
savefig(p0,"1_0.png")

p1 = VMD.plot(v,k=1)
savefig(p1,"1_1.png")

p2 = VMD.plot(v,k=2)
savefig(p2,"1_2.png")

p3 = VMD.plot(v,k=3)
savefig(p3,"1_3.png")
```


## plot the original signal and spectrum
![](1_0.png)

## plot the 1st decomposed signal and spectrum
![](1_1.png)

## plot the 2st signal and spectrum
![](1_2.png)

## plot the 3st decomposed signal and spectrum
![](1_3.png)
