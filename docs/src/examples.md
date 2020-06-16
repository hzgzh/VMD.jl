# Examples

```@repl
using VMD,Random,Plots
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
K = 4;              # 3 modes
DC = false;             # no DC part imposed
init = 1;           # initialize omegas uniformly
tol = 1e-7;


v = vmd(f ; alpha = alpha,tau = tau,K = K,DC = false,init = 1,tol = tol,sample_frequency = sample_frequency)

# the first mode frequency
print("1st mode frequency $(n_mode(v,1))")

p0 = VMD.plot(v,k=0)
savefig(p0,"0.png")

p1 = VMD.plot(v,k=1)
savefig(p1,"1.png")

p2 = VMD.plot(v,k=2)
savefig(p2,"2.png")

p3 = VMD.plot(v,k=3)
savefig(p3,"3.png")

```

## plot the original signal and spectrum
![](0.png)

## plot the 1st decomposed signal and spectrum
![](1.png)

## plot the 2st signal and spectrum
![](2.png)

## plot the 3st decomposed signal and spectrum
![](3.png)

