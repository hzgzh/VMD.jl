var documenterSearchIndex = {"docs":
[{"location":"api/#Reference-guide-1","page":"Reference","title":"Reference guide","text":"","category":"section"},{"location":"api/#Index-1","page":"Reference","title":"Index","text":"","category":"section"},{"location":"api/#","page":"Reference","title":"Reference","text":"","category":"page"},{"location":"api/#Public-API-1","page":"Reference","title":"Public API","text":"","category":"section"},{"location":"api/#","page":"Reference","title":"Reference","text":"Modules = [VMD]\nOrder   = [:type,:function]","category":"page"},{"location":"api/#VMD.Vmd","page":"Reference","title":"VMD.Vmd","text":"Vmd{T,S<:Int}\n\nVariational Mode Decompostion \nReturn by vmd function.\n\n** Arguments **\n\nsignal::Array{T,1}:signal\nK::S              :modes \nsignal_d::Array{T,Val{K}}:decomposed signals\nfreqs::Array{T,Val{K}}:decomposed spectrums\nsamples::S:signal sampled frequency\n\n\n\n\n\n","category":"type"},{"location":"api/#VMD.compare-Tuple{Vmd}","page":"Reference","title":"VMD.compare","text":"compare(v::Vmd)\n\nvisual compare origin signal and sum of decomposed signal\n\n\n\n\n\n","category":"method"},{"location":"api/#VMD.n_component-Tuple{Vmd,Any}","page":"Reference","title":"VMD.n_component","text":"n_component(v::Vmd,k)\n\nreturn No k component decomposed signal\n\n\n\n\n\n","category":"method"},{"location":"api/#VMD.n_mode-Tuple{Vmd,Any}","page":"Reference","title":"VMD.n_mode","text":"n_mode(v::Vmd,k)\nreturn the k mode frequency\n\n\n\n\n\n","category":"method"},{"location":"api/#VMD.plot-Tuple{Vmd}","page":"Reference","title":"VMD.plot","text":"plot(vmd::Vmd;k=1)\n\nvisual the decomposed signals and spectrums\n\n** Argument **\n\nvmd::Vmd : vmd\nk::Int   : 0-original signal 1-first component enforcement\n\n\n\n\n\n","category":"method"},{"location":"api/#VMD.vmd-Union{Tuple{Array{Typ,1}}, Tuple{Typ}} where Typ","page":"Reference","title":"VMD.vmd","text":"vmd(signal;alpha=100, tau=0, K=3, DC=false, init=1, tol=1e-6,samples=100)\n\nVariational Mode Decomposition\nReturn Vmd\n\n** Argument **\n\nsignal: the time domain signal (1D) to be decomposed\nalpha : the balancing parameter of the data-fidelity constraint\ntau   : time-step of the dual ascent ( pick 0 for noise-slack )\nK     : the number of modes to be recovered\nDC    : true if the first mode is put and kept at DC (0-freq)\ninit  : 0 = all omegas start at 0           1 = all omegas start uniformly distributed           2 = all omegas initialized randomly\ntol   : tolerance of convergence criterion; typically around 1e-6\nsample_frequency : samples frequency(eg:100/s)\n\n** Example **\n\nT = 1000;\nt = (1:T)/T;\nsample_frequency = 1000;\n\n# center frequencies of components\nf_1 = 2;\nf_2 = 24;\nf_3 = 288;\n\n# modes\nv_1 = @. cos(2*pi*f_1*t);\nv_2 = @. 1/4*(cos(2*pi*f_2*t));\nv_3 = @. 1/16*(cos(2*pi*f_3*t));\n\n# composite signal, including noise\nf = v_1 + v_2 + v_3 + 0.1*randn(length(v_1));\n\n# some sample parameters for VMD\nalpha = 2000;       # moderate bandwidth constraint\ntau = 0;            # noise-tolerance (no strict fidelity enforcement)\nK = 3;              # 3 modes\nDC = false;             # no DC part imposed\ninit = 1;           # initialize omegas uniformly\ntol = 1e-7;\n\n\nv = vmd(f ; alpha = alpha,tau = tau,K = K,DC = false,init = 1,tol = tol,sample_frequency = sample_frequency)\n\n# plot original signal and spectrum\nplot(v;k = 0)\n\n# plot first decomposed component and spectrum\nplot(v;k = 1)\n\n\n\n\n\n","category":"method"},{"location":"example2/#Example2-1","page":"-","title":"Example2","text":"","category":"section"},{"location":"example2/#","page":"-","title":"-","text":"The second example uses a quadratic trend, a chirp signal, and a third mode with sharp transition between two constant frequencies","category":"page"},{"location":"example2/#","page":"-","title":"-","text":"$","category":"page"},{"location":"example2/#","page":"-","title":"-","text":"f_{\\mathrm{Sig} 2}(t)=6 t^{2}+\\cos \\left(10 \\pi t+10 \\pi t^{2}\\right)+\\left{\\begin{array}{ll} \\cos (60 \\pi t) & t \\leq 0.5 \\\n\\cos (80 \\pi t-10 \\pi) & t>0.5 \\end{array}\\right. $","category":"page"},{"location":"example2/#","page":"-","title":"-","text":"For t in 01the chirp’s instantaneous frequency varies linearly between 10pi an 30pi.Consequently, the theoretical center frequency of the mode is located at 20pi . The piecewise-constant bi-harmonic has spectral peaks expected at 60pi and 80pi","category":"page"},{"location":"example2/#","page":"-","title":"-","text":"Again, with VMD the estimated center frequencies con- verge to the expected frequencies precisely. Here, we chose to decompose into four modes, thus assigning each half of the piecewise-constant frequency signal to a separate mode. The spectral partitioning can be nicely appreciated in the spectral plot of the different modes. Here, EMD does a better job recov- ering the quadratic trend, and, correspondingly, the first oscilla- tion. However, EMD is unable to separate the two pieces of the piecewise constant frequency signal. .","category":"page"},{"location":"example2/#","page":"-","title":"-","text":"using VMD,Random,Plots\nT = 1000;\nt = (1:T)/T;\nsample_frequency = 1000;\n\n# center frequencies of components\nf_1 = 10;\nf_2 = 60;\nf_3 = 80;\n\n# modes\nv_1 = @. 6t^2;\nv_2 = @. cos(10π*t+10π*t^2)\nv_3 = [t1>=0.5 ? cos(60π*t1) : cos(80π*t1-10π) for t1 in t] \n\n# composite signal, including noise\nf = v_1 + v_2 + v_3 + 0.1*randn(length(v_1));\n\n# some sample parameters for VMD\nalpha = 2000;       # moderate bandwidth constraint\ntau = 0;            # noise-tolerance (no strict fidelity enforcement)\nK = 4;              # 3 modes\nDC = false;             # no DC part imposed\ninit = 0;           # initialize omegas uniformly\ntol = 1e-7;\n\n\nv = vmd(f ; alpha = alpha,tau = tau,K = K,DC = false,init = 1,tol = tol,sample_frequency = sample_frequency)\n\n# the first mode frequency\nprint(\"1st mode frequency $(n_mode(v,1))\")\n\np0 = VMD.plot(v,k=0)\nsavefig(p0,\"2_0.png\")\n\np1 = VMD.plot(v,k=1)\nsavefig(p1,\"2_1.png\")\n\np2 = VMD.plot(v,k=2)\nsavefig(p2,\"2_2.png\")\n\np3 = VMD.plot(v,k=3)\nsavefig(p3,\"2_3.png\")\n","category":"page"},{"location":"example2/#plot-the-original-signal-and-spectrum-1","page":"-","title":"plot the original signal and spectrum","text":"","category":"section"},{"location":"example2/#","page":"-","title":"-","text":"(Image: )","category":"page"},{"location":"example2/#plot-the-1st-decomposed-signal-and-spectrum-1","page":"-","title":"plot the 1st decomposed signal and spectrum","text":"","category":"section"},{"location":"example2/#","page":"-","title":"-","text":"(Image: )","category":"page"},{"location":"example2/#plot-the-2st-signal-and-spectrum-1","page":"-","title":"plot the 2st signal and spectrum","text":"","category":"section"},{"location":"example2/#","page":"-","title":"-","text":"(Image: )","category":"page"},{"location":"example2/#plot-the-3st-decomposed-signal-and-spectrum-1","page":"-","title":"plot the 3st decomposed signal and spectrum","text":"","category":"section"},{"location":"example2/#","page":"-","title":"-","text":"(Image: )","category":"page"},{"location":"example3/#Example-3-1","page":"-","title":"Example 3","text":"","category":"section"},{"location":"example3/#","page":"-","title":"-","text":"The third synthetic signal has intrawave fre- quency modulation:","category":"page"},{"location":"example3/#","page":"-","title":"-","text":"$","category":"page"},{"location":"example3/#","page":"-","title":"-","text":"f_{\\mathrm{Sig} 3}(t)=\\frac{1}{1.2+\\cos (2 \\pi t)}+\\frac{\\cos (32 \\pi t+0.2 \\cos (64 \\pi t))}{1.5+\\sin (2 \\pi t)} $","category":"page"},{"location":"example3/#","page":"-","title":"-","text":"The signal, its three constituent modes, and the composite Fourier spectrum are shown in Fig. 11. While the first, bell-shaped component has mostly low-pass content, the second mode’s main peak is clearly identified at 32pi . However, due to the non-linear intrawave frequency modulation, an important amount of higher-order harmonics are also observed. This second component obviously violates the narrow-band assumption, and one would naturally expect some difficulties recovering this mode using VMD.","category":"page"},{"location":"example3/#","page":"-","title":"-","text":"The corresponding VMD results are illustrated in Fig. 11. The non-zero omega_2 quickly converges to the correct main frequency 32pi . The higher order harmonics are not uniquely attributed to the second mode, however, but shared between both modes. Consequently, the intrawave frequency modulation is shared by both modes, creating some ripples in the otherwise low-fre- quency mode.","category":"page"},{"location":"example3/#","page":"-","title":"-","text":"using VMD,Random,Plots\nT = 1000;\nt = (1:T)/T;\nsample_frequency = 1000;\n\n# modes\nv_1 = @. 1.0/(1.2+cos(2π*t))\nv_2 = @. cos(32π*t+0.2cos(64π*t))\nv_3 = @. 1.5+sin(2π*t)\n\n# composite signal, including noise\nf = v_1 + v_2./v_3 + 0.1*randn(length(v_1));\n\n# some sample parameters for VMD\nalpha = 2000;       # moderate bandwidth constraint\ntau = 0;            # noise-tolerance (no strict fidelity enforcement)\nK = 4;              # 3 modes\nDC = false;             # no DC part imposed\ninit = 0;           # initialize omegas uniformly\ntol = 1e-7;\n\n\nv = vmd(f ; alpha = alpha,tau = tau,K = K,DC = false,init = init,tol = tol,sample_frequency = sample_frequency)\n\n# the first mode frequency\nprint(\"1st mode frequency $(n_mode(v,1))\")\n\np0 = VMD.plot(v,k=0)\nsavefig(p0,\"3_0.png\")\n\np1 = VMD.plot(v,k=1)\nsavefig(p1,\"3_1.png\")\n\np2 = VMD.plot(v,k=2)\nsavefig(p2,\"3_2.png\")\n\np3 = VMD.plot(v,k=3)\nsavefig(p3,\"3_3.png\")","category":"page"},{"location":"example3/#plot-the-original-signal-and-spectrum-1","page":"-","title":"plot the original signal and spectrum","text":"","category":"section"},{"location":"example3/#","page":"-","title":"-","text":"(Image: )","category":"page"},{"location":"example3/#plot-the-1st-decomposed-signal-and-spectrum-1","page":"-","title":"plot the 1st decomposed signal and spectrum","text":"","category":"section"},{"location":"example3/#","page":"-","title":"-","text":"(Image: )","category":"page"},{"location":"example3/#plot-the-2st-signal-and-spectrum-1","page":"-","title":"plot the 2st signal and spectrum","text":"","category":"section"},{"location":"example3/#","page":"-","title":"-","text":"(Image: )","category":"page"},{"location":"example3/#plot-the-3st-decomposed-signal-and-spectrum-1","page":"-","title":"plot the 3st decomposed signal and spectrum","text":"","category":"section"},{"location":"example3/#","page":"-","title":"-","text":"(Image: )","category":"page"},{"location":"examples/#Examples-1","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/#","page":"Examples","title":"Examples","text":"using VMD,Random,Plots\nT = 1000;\nt = (1:T)/T;\nsample_frequency = 1000;\n\n# center frequencies of components\nf_1 = 2;\nf_2 = 24;\nf_3 = 288;\n\n# modes\nv_1 = @. cos(2*pi*f_1*t);\nv_2 = @. 1/4*(cos(2*pi*f_2*t));\nv_3 = @. 1/16*(cos(2*pi*f_3*t));\n\n# composite signal, including noise\nf = v_1 + v_2 + v_3 + 0.1*randn(length(v_1));\n\n# some sample parameters for VMD\nalpha = 2000;       # moderate bandwidth constraint\ntau = 0;            # noise-tolerance (no strict fidelity enforcement)\nK = 3;              # 3 modes\nDC = false;             # no DC part imposed\ninit = 1;           # initialize omegas uniformly\ntol = 1e-7;\n\n\nv = vmd(f ; alpha = alpha,tau = tau,K = K,DC = false,init = 1,tol = tol,sample_frequency = sample_frequency)\n\n# the first mode frequency\nprint(\"1st mode frequency $(n_mode(v,1))\")\n\np0 = VMD.plot(v,k=0)\nsavefig(p0,\"0.png\")\n\np1 = VMD.plot(v,k=1)\nsavefig(p1,\"1.png\")\n\np2 = VMD.plot(v,k=2)\nsavefig(p2,\"2.png\")\n\np3 = VMD.plot(v,k=3)\nsavefig(p3,\"3.png\")\n","category":"page"},{"location":"examples/#plot-the-original-signal-and-spectrum-1","page":"Examples","title":"plot the original signal and spectrum","text":"","category":"section"},{"location":"examples/#","page":"Examples","title":"Examples","text":"(Image: )","category":"page"},{"location":"examples/#plot-the-1st-decomposed-signal-and-spectrum-1","page":"Examples","title":"plot the 1st decomposed signal and spectrum","text":"","category":"section"},{"location":"examples/#","page":"Examples","title":"Examples","text":"(Image: )","category":"page"},{"location":"examples/#plot-the-2st-signal-and-spectrum-1","page":"Examples","title":"plot the 2st signal and spectrum","text":"","category":"section"},{"location":"examples/#","page":"Examples","title":"Examples","text":"(Image: )","category":"page"},{"location":"examples/#plot-the-3st-decomposed-signal-and-spectrum-1","page":"Examples","title":"plot the 3st decomposed signal and spectrum","text":"","category":"section"},{"location":"examples/#","page":"Examples","title":"Examples","text":"(Image: )","category":"page"},{"location":"example1/#Example1-1","page":"-","title":"Example1","text":"","category":"section"},{"location":"example1/#","page":"-","title":"-","text":"The first signal is a composition of three simple components, namely a general linear trend and two dif- ferent harmonics: f_mathrmSig 1(t)=6 t+cos (8 pi t)+05 cos (40 pi t)","category":"page"},{"location":"example1/#","page":"-","title":"-","text":"The signal, its three constituent modes, and the composite Fourier spectrum are shown in Fig. 9. The main challenge of this signal is the linear growth term, whose higher order harmonics spread over the whole spectrum.","category":"page"},{"location":"example1/#","page":"-","title":"-","text":"The recovered VMD modes constitute a nice partition of the input spectrum, with each mode being clearly dominant around its respective center frequency. The three modes in time domain show nice separation into three distinct signals of characteristic oscillations. Compared to the EMD results, the proposed VMD method performs clearly better, with less spurious oscillations in the trend and mid-frequency signal.","category":"page"},{"location":"example1/#","page":"-","title":"-","text":"using VMD,Random,Plots\nT = 1000;\nt = (1:T)/T;\nsample_frequency = 1000;\n\n# modes\nv_1 = @. 6t;\nv_2 = @. cos(8π*t)\nv_3 = @. 0.5cos(40π*t) \n\n# composite signal, including noise\nf = v_1 + v_2 + v_3 + 0.1*randn(length(v_1));\n\n# some sample parameters for VMD\nalpha = 2000;       # moderate bandwidth constraint\ntau = 0;            # noise-tolerance (no strict fidelity enforcement)\nK = 3;              # 3 modes\nDC = false;             # no DC part imposed\ninit = 0;           # initialize omegas uniformly\ntol = 1e-7;\n\n\nv = vmd(f ; alpha = alpha,tau = tau,K = K,DC = false,init = init,tol = tol,sample_frequency = sample_frequency)\n\n\np0 = VMD.plot(v,k=0)\nsavefig(p0,\"1_0.png\")\n\np1 = VMD.plot(v,k=1)\nsavefig(p1,\"1_1.png\")\n\np2 = VMD.plot(v,k=2)\nsavefig(p2,\"1_2.png\")\n\np3 = VMD.plot(v,k=3)\nsavefig(p3,\"1_3.png\")","category":"page"},{"location":"example1/#plot-the-original-signal-and-spectrum-1","page":"-","title":"plot the original signal and spectrum","text":"","category":"section"},{"location":"example1/#","page":"-","title":"-","text":"(Image: )","category":"page"},{"location":"example1/#plot-the-1st-decomposed-signal-and-spectrum-1","page":"-","title":"plot the 1st decomposed signal and spectrum","text":"","category":"section"},{"location":"example1/#","page":"-","title":"-","text":"(Image: )","category":"page"},{"location":"example1/#plot-the-2st-signal-and-spectrum-1","page":"-","title":"plot the 2st signal and spectrum","text":"","category":"section"},{"location":"example1/#","page":"-","title":"-","text":"(Image: )","category":"page"},{"location":"example1/#plot-the-3st-decomposed-signal-and-spectrum-1","page":"-","title":"plot the 3st decomposed signal and spectrum","text":"","category":"section"},{"location":"example1/#","page":"-","title":"-","text":"(Image: )","category":"page"},{"location":"#VMD.jl-Document-1","page":"Home","title":"VMD.jl Document","text":"","category":"section"},{"location":"#Introduction-1","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"The goal of VMD is to decompose a real valued input signal f into a discrete number of sub-signals (modes),u_k , that have specific sparsity properties while reproducing the input Here, the sparsity prior of each mode is chosen to be its bandwidth in spectral domain. In other words, we assume each mode k to be mostly compact around a center pulsation omega_k , which is to be determined along with the decomposition.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"CurrentModule = VMD","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Pages   = [\"index.md\",\"examples.md\"]\nModules = [VMD]\nOrder   = [:function, :type]","category":"page"}]
}