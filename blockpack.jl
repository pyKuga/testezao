function F2DOF(dt,Param)

    m1 = Param[1];
    m2 = Param[2];

    c1 = Param[3];
    c2 = Param[4];
    c3 = Param[5];

    k1 = Param[6];
    k2 = Param[7];
    k3 = Param[8];

    M = [m1 0; 0 m2];
    C = [c1+c2 -c2; -c2 c2+c3];
    K = [k1+k2 -k2; -k2 k2+k3];

    Id = [1 0; 0 1];
    Z = [0 0; 0 0];

    invM = inv(M);

    A = [Z Id; -invM*K -invM*C];
    B = [Z; invM];
    B = B[:,1];
    C = [1 0 0 0; 0 1 0 0];
    #D = [0 0;0 0];

    sistema = ss(A,B,C,0);
    sistemaD = c2d(sistema,dt);

    return sistema, sistemaD 
end  
    
#gera um sinal normal normatizado manualmente
function noise_gen(seed,r,c)
    merst = MersenneTwister(seed); #mersenne twister generator
    signal = randn(merst, Float64,(r,c)); #gaussian noise
    mu = mean(signal,dims=2);
    signal_mu = signal.-mu;
    across = (maximum(signal.-mu,dims=2) - minimum(signal.-mu,dims=2))/2;
    return signal_mu./across
end