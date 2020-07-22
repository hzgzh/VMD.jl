## Example2

The second example uses a quadratic trend,
a chirp signal, and a third mode with sharp transition between
two constant frequencies

$$
f_{\mathrm{Sig} 2}(t)=6 t^{2}+\cos \left(10 \pi t+10 \pi t^{2}\right)+\left\{\begin{array}{ll}
\cos (60 \pi t) & t \leq 0.5 \\
\cos (80 \pi t-10 \pi) & t>0.5
\end{array}\right.
$$

For $t \in [0,1]$the chirp’s instantaneous frequency varies linearly between $10\pi$ an $30\pi$.Consequently, the theoretical center frequency of the mode is located at $20\pi$
. The piecewise-constant bi-harmonic has
spectral peaks expected at $60\pi$ and $80\pi$

Again, with VMD the estimated center frequencies
con-
verge to the expected frequencies precisely. Here, we chose to
decompose into four modes, thus assigning each half of the
piecewise-constant frequency signal to a separate mode. The
spectral partitioning can be nicely appreciated in the spectral
plot of the different modes. Here, EMD does a better job recov-
ering the quadratic trend, and, correspondingly, the first oscilla-
tion. However, EMD is unable to separate the two pieces of the
piecewise constant frequency signal.
.
```@repl
using VMD,Random,Plots
T = 1000;
t = (1:T)/T;
sample_frequency = 1000;

# center frequencies of components
f_1 = 10;
f_2 = 60;
f_3 = 80;

# modes
v_1 = @. 6t^2;
v_2 = @. cos(10π*t+10π*t^2)
v_3 = [t1>=0.5 ? cos(60π*t1) : cos(80π*t1-10π) for t1 in t] 

# composite signal, including noise
f = v_1 + v_2 + v_3 + 0.1*randn(length(v_1));

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
savefig(p0,"2_0.png")

p1 = VMD.plot(v,k=1)
savefig(p1,"2_1.png")

p2 = VMD.plot(v,k=2)
savefig(p2,"2_2.png")

p3 = VMD.plot(v,k=3)
savefig(p3,"2_3.png")

p4 = VMD.plot(v,k=4)
savefig(p4,"2_4.png")
```
## plot the original signal and spectrum
![](2_0.png)

## plot the 1st decomposed signal and spectrum
![](2_1.png)

## plot the 2st signal and spectrum
![](2_2.png)

## plot the 3st decomposed signal and spectrum
![](2_3.png)

## plot the 4st decomposed signal and spectrum
![](2_4.png)