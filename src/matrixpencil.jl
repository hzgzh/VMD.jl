using LinearAlgebra
using Random

struct MatrixPencil{T}
    y::Array{T,1}
    freq::Array{T,1}
    dr::Array{T,1}
    phase::Array{T,1}
    amplitude::Array{T,1}
    sample_frequency::T
    function MatrixPencil(y,freq,dr,phase,amplitude,sample_frequency)
        new{eltype(y)}(y,freq,dr,phase,amplitude,sample_frequency)
    end
end

function Base.show(io::IO,mp::MatrixPencil{T}) where T
    println("frequencies  (Hz) : $(mp.freq)")
    println("damp ratios  (%)  : $(100mp.dr)")
    println("phases       (°)  : $(mp.phase)")
    println("amplitudes        : $(mp.amplitude)")
end

frequencies(mp::MatrixPencil) = mp.freq
dampratio(mp::MatrixPencil) = mp.dr
phase(mp::MatrixPencil) = mp.phase
amplitude(mp::MatrixPencil) = mp.amplitude

"""
    matrixpencil(y,yz,dt)
    矩阵束方法分析低频振荡

- `y`        ：     输入样本数据,必须是列向量
- `threshold`：     阶数估计阀值
- `dt`       ：     时间间隔
"""
function matrixpencil(y,threshold,dt)
    N=length(y)#样本长度
    L=round(Int,N*0.3)#矩阵束参数
    
    bz=1

    #第一步：构造Hankel矩阵
    Y=zeros(N-L,L+1)
    for i=1:L+1,j=1:N-L
        Y[j,i]=y[i+j-1]
    end

    #第二步：奇异值分解确定模态阶
    U,D,V=svd!(Y)
    
    i=0
    while bz>threshold
        i=i+1
        bz=D[i]/D[1]
    end
    M=i-1

    #第三步：构造D'矩阵
    D =diagm(D)
    D1=view(D,:,1:M)

    #第四步：求取模态的阻尼和频率
    V1=view(V,1:L,1:M) #截取V前M个主导右特征向量的第1行到第L行
    V2=view(V,2:L+1,1:M)#截取V前M个主导右特征向量的第2行到第L+1行
    Y1=U*D1*V1'
    Y2=U*D1*V2'
    #Y1w为Y1伪逆矩阵%% Y1w=inv(Y1'*Y1)*A';
    Y1w=pinv(Y1)#(Y1'*Y1)\Y1';%Y1w为Y1伪逆矩阵
    G=Y1w*Y2
    
    #求取G的非零特征值，再求阻尼和频率
    fvalue = eigvals!(G)
    fvalue = reverse(fvalue)
    
    f=zeros(M)
    da=zeros(M)
    c=zeros(Complex{Float64},M)
    for i=1:M
        if abs(fvalue[i])<1e-4
             da[i]=0  
             f[i]=0
        else
            println(fvalue[i])
            c[i]=log(fvalue[i])/dt
            println(c[i])
            da[i]=real(c[i])  
            f[i]=imag(c[i])/(2π)
        end
    end

    #振荡幅值的求解
    z=fvalue[1:M]
    amplitude=zeros(M)
    phase=zeros(M)
    #范德蒙德矩阵
    Z=zeros(Complex{Float64},N,M)
    for i=1:N,j=1:M
            Z[i,j]=z[j]^(i-1)
    end
    #范德蒙德矩阵求逆
    Zw=pinv(Z)
    #
    b=Zw*y
    for i=1:M
        amplitude[i]=abs(b[i])#第一列为各次幅值
        phase[i]=atan(imag(b[i])/real(b[i]))*180/pi #第二列为各次初相位
        #初始相位校正
        if real(b[i])<0
            phase[i]=180+phase[i]
        end
    end
    fv = [f.>0]
    MatrixPencil(y,f[fv...],da[fv...],phase[fv...],2amplitude[fv...],1/dt)
end

using Plots

function sim(mp::MatrixPencil)
    n = length(mp.y)
    t = 0:1/mp.sample_frequency:1/mp.sample_frequency*(n-1)
    f = frequencies(mp)
    dr = dampratio(mp)
    ps = phase(mp)
    am = amplitude(mp)
    ys = zeros(length(mp.y))
    for i in 1:length(f)
        ys += am[i]*exp.(-abs(dr[i])*t).*cos.(2π*f[i]*t.+ps[1]*π/180)
    end
    plot(t,y,xlabel="time s",label="actual signal")
    plot!(t,ys,label="sim signal")
end

    # Write your tests here.
    dt =0.01
    t = 0:dt:dt*1000
    sample_frequency = 100

# center frequencies of components
    f_1 = 2
    f_2 = 24
    f_3 = 48
    

# modes
    v_1 = @. 10exp(-0.2*t)cos(2 * pi * f_1 * t+π/2)
    v_2 = @. 4exp(-0.1t) * (cos(2 * pi * f_2 * t+π/4))
    v_3 = @. 2exp(-0.05t) * (cos(2 * pi * f_3 * t))

# composite signal, including noise
    y = v_1 + v_2 + v_3 + 0.0001 * randn(length(v_1))
    #y = 0.01randn(1000)
    dt = 1/sample_frequency

    @time mp = matrixpencil(y,0.8,dt)