# VMD.jl Document


## Introduction
The goal of VMD is to decompose a real valued input signal $f$
into a discrete number of sub-signals (modes),$u_k$
, that have
specific sparsity properties while reproducing the input Here,
the sparsity prior of each mode is chosen to be its bandwidth
in spectral domain. In other words, we assume each mode $k$
to
be mostly compact around a center pulsation $\omega_k$
, which is to be
determined along with the decomposition.

```@meta
CurrentModule = VMD
```

```@index
Pages   = ["index.md","examples.md"]
Modules = [VMD]
Order   = [:function, :type]
```



