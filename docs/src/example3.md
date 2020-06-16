### Example 3

The third synthetic signal has intrawave fre-
quency modulation:

$$
f_{\mathrm{Sig} 3}(t)=\frac{1}{1.2+\cos (2 \pi t)}+\frac{\cos (32 \pi t+0.2 \cos (64 \pi t))}{1.5+\sin (2 \pi t)}
$$

The signal, its three constituent modes, and the composite
Fourier spectrum are shown in Fig. 11. While the first,
bell-shaped component has mostly low-pass content, the
second mode’s main peak is clearly identified at $32\pi$
. However,
due to the non-linear intrawave frequency modulation, an
important amount of higher-order harmonics are also observed.
This second component obviously violates the narrow-band
assumption, and one would naturally expect some difficulties
recovering this mode using VMD.

The corresponding VMD results are illustrated in Fig. 11. The
non-zero $\omega_2$
quickly converges to the correct main frequency $32\pi$
. The higher order harmonics are not uniquely attributed
to the second mode, however, but shared between both modes.
Consequently, the intrawave frequency modulation is shared by
both modes, creating some ripples in the otherwise low-fre-
quency mode.

```@repl
using VMD,Random,Plots
T = 1000;
t = (1:T)/T;
sample_frequency = 1000;

# modes
v_1 = @. 1.0/(1.2+cos(2π*t))
v_2 = @. cos(32π*t+0.2cos(64π*t))
v_3 = @. 1.5+sin(2π*t)

# composite signal, including noise
f = v_1 + v_2./v_3 + 0.1*randn(length(v_1));

# some sample parameters for VMD
alpha = 2000;       # moderate bandwidth constraint
tau = 0;            # noise-tolerance (no strict fidelity enforcement)
K = 4;              # 3 modes
DC = false;             # no DC part imposed
init = 0;           # initialize omegas uniformly
tol = 1e-7;


v = vmd(f ; alpha = alpha,tau = tau,K = K,DC = false,init = init,tol = tol,sample_frequency = sample_frequency)

# the first mode frequency
print("1st mode frequency $(n_mode(v,1))")

p0 = VMD.plot(v,k=0)
savefig(p0,"3_0.png")

p1 = VMD.plot(v,k=1)
savefig(p1,"3_1.png")

p2 = VMD.plot(v,k=2)
savefig(p2,"3_2.png")

p3 = VMD.plot(v,k=3)
savefig(p3,"3_3.png")
```

## plot the original signal and spectrum
![](3_0.png)

## plot the 1st decomposed signal and spectrum
![](3_1.png)

## plot the 2st signal and spectrum
![](3_2.png)

## plot the 3st decomposed signal and spectrum
![](3_3.png)