using Flux
using CUDA


#School of Mechanical Engineering - State University of Campinas
#Paulo Yoshio Kuga

#This package is beign written for "Análise de Robustez em Métodos de Identificação de Sistemas", a undergraduate ressearch for PIBIC.
#This package is also being used in "Análise de Robustez em Métodos de Identificação de Sistemas e Redes Neurais" the second signed PIBIC propposal.

#Main content: sysid algorithms and it's support functions

#------------------------------------------GENERAL FUNCTIONS-------------------------------------------------

function CheckData(Y,ny,ns)
    if ns < ny
        Y = Y';
        ny,ns = size(Y);
    end
    return Y, ny, ns
end

#------------------------------------------ERA RELATED FUNCTIONS-------------------------------------------------

function Hankel(h,p,ny)
    h = h[:,2:p+1];
    ch = round(Int,p/2); #colunas da matriz de hankel
    rh = ch; #linhas da matriz de hankel
    H = zeros(rh*ny,ch+1); #estabelecer matriz de hankel baseada no número de linhas e entradas desejadas
    for i in 1:rh #preenche a matriz de hankel (dificil explicar a logica, mas em suma tenta tomar os espaços já pré-determinados para preencher a matriz)
        H[(1:ny).+((i-1)*ny),:] = h[:,(1:ch+1).+(i-1)]
    end
    H1 = H[:,1:ch];
    H2 = H[:,(1:ch).+1];
    return H1, H2
end

function ERA_K(Y,nx,nu,dt,p,ny)
    H1,H2 = Hankel(Y,p,ny);
    U, S, V = svd(H1); #faz a SVD
    S = diagm(S[1:nx]); #transforma os autovalores do svd em uma matriz diagonal
    U = U[:,1:nx]; #toma as colunas relativas a ordem do modelo
    V = V[:,1:nx]; #confirmar com o prof a hipótese do V ser tomado como coluna devido a transposição e da definição do SVD em si
    G = S.^0.5; #tira a raiz quadrada dos autovalores
    Q = inv(G); #inverte-os
    Or = U*G; #observabilidade do modelo
    Cr = G*V'; #controlabilidade do modelo
    A = Q*U'*H2*V*Q; #gera a matriz A
    C = Or[1:ny,:]; #com base no número de entradas do modelo, gera o 
    B = Cr[:,1:nu]; #numero de entradas é um dado do problema
    D = Y[:,1];
    return ss(A,B,C,D,dt)
end

#Y: the system impulse response 
#nx: model order
#nu: number of inputs
#dt: Y sample period
function ImpulseERA(Y,nx,nu,dt,p)
    ny,ns = size(Y);
    Y,ny,ns = CheckData(Y,ny,ns);
    return ERA_K(Y,nx,nu,dt,p,ny)
end

#bulid the U matrix presented in Dirceu's Thesis.
function BulidU(U,p,nu,ns) 
    U = U[:,1:ns-1];
    cut = ns-p  #cut parameter for lines
    Hu = zeros(nu*cut,ns-1);
    for i=1:cut
        Hu[(1:nu).+((i-1)*nu),i:ns-1] = U[:,1:ns-1-(i-1)];
    end
    return Hu
end

function GRA(Y,U,nx,dt,p)
    ny,ns = size(Y);
    nu, = size(U); #yes, for now, i'm putting some faith in user.
    Y,ny,ns = CheckData(Y,ny,ns);
    U,nu,ns = CheckData(U,nu,ns);

    Uh = BulidU(U,p,nu,ns);
    Yh = Y[:,1:ns-1];

    h = Yh*Uh'*inv(Uh*Uh');
    return ERA_K(h,nx,nu,dt,p,ny)
    

end


#------------------------------------------ARX RELATED FUNCTIONS-------------------------------------------------

#Hankel pré-determinado
#np ordem do polinimio
#ns numero de amostras
function Hankel_PD(Y,np,linhas)
    H = zeros(linhas,np); #estabelecer matriz de hankel baseada no número de linhas e entradas desejadas
    ny,ns = size(Y);
    Yf = reshape(Y,(ny*ns,1));
    for i in 1:np
        H[:,i] = Yf[(1:linhas).+(i-1)*ny]#reshape(Y[:,(1:(ns-np)).+(i-1)],linhas,1);
    end
    return H
end

function Hankel_PD2(U,np,lines)
    nu,ns = size(U);
    Ut = reshape(U,(1,ns*nu));
    F = zeros(lines,np*nu); #estabelecer matriz de hankel baseada no número de linhas e entradas desejadas
    for i in 1:lines 
        F[i,:] = Ut[(1:np*nu).+(i-1)*nu];
    end
    return F
end

#Y: multiple output signal
#U: multiple input signal
#na: denominator polynomial order na>=nb
#nb: numerator polynomial order nb>=1
function ARX_K(Y,U, na, nb,dt)
    ny,ns = size(Y); #determines the number of outputs and the number of samples
    nu,_ = size(U); #determines the number of inputs
    #if the signal or the input come as rows being the samples, it changes the order for columns being samples
    Y,ny,ns = CheckData(Y,ny,ns); 
    U,nu,ns = CheckData(U,nu,ns);

    linhas = (ns-na)*ny; #since we want na coefficients to our denominator, we need to organize the samples.
    colunas = na+nb*ny*nu; #it's related to the number of coefficients. 
    H = zeros(linhas,colunas); #hankel matrix is reserved in memory

    #for building our StateSpace system, we storage a matrix with the TF discrete type, to populate the array
    G = Matrix{TransferFunction{Discrete{Float64}, ControlSystemsBase.SisoRational{Float64}}}(undef,ny,nu); 
    

    #professor kurka formulation for MIMO-ARX with least square
    Id = Matrix{Float64}(I,ny,ny);
    H[:,1:nb*ny*nu] = kron(Hankel_PD2(U,nb,ns-na),Id);
    H[:,(1:na).+nb*ny*nu] = -Hankel_PD(Y,na,linhas);
    b = reshape(Y[:,na+1:ns],linhas,1);

    Coef = H'*H\(H'*b);
    A = reverse(push!(Coef[(1:na).+nb*ny],1));
    tfA = tf(1,A,dt);


    B = Coef[1:(nb*ny)];
 

    #TF functions assembly and asssignment to the G matrix (the TF function one)
    Bc = zeros(ny,nb+1); #a new matrix is called 
    Bc[:,1:nb] = reverse(reshape(B,(ny*nu,nb)),dims=2); #it gets the coefficients and separe them into a nb row matrix, nu*ny columns. then the matrix is transposed, and reversed
    for i=1:ny
        for j = 1:nu
            G[i,j] = 1/tf(1,Bc[i+j-1,:],dt)*tfA;
        end        
    end

    

    return array2mimo(G); #conversion to a mimo system
    #returns the space state form of the system, not the TF.
end


#------------------------------------------NARX RELATED FUNCTIONS-------------------------------------------------

# mutable struct NNARXModel
#     nn::Chain
#     na::Int64
#     nb::Int64
#     nu::Int64
#     ny::Int64
# end

# function Factory(depth,large)
#     return [ny*na+nu*nb; rand(1:large,depth);2]
# end


# function NARX(ny,nu,na,nb,V)
#     nv, = size(V)
#     F = [[tanh for i=1:nv-2]; identity] |> gpu
#     neurons = [Dense(V[i],V[i+1],F[i]) for i=1:nv-1] |> gpu
#     model = Chain(neurons) |> gpu

#     return NNARXModel(model,na,nb,nu,ny)
# end


# function DataNARX(Y,U, ny, nu, na, nb,ns)
    
#     Yt = [vec(reshape(Y[:,i],(ny,1))) for i = 2:ns] |> gpu
    
#     r = ns - na
#     Yr = DataNorm([zeros(ny,na-1) Y]); #we consider that the system for t < 0 is at rest   
#     U = DataNorm([zeros(nu,nb-1) U]);
#     H = zeros(na*ny+nu*nb,ns-1); #this is the matrix who allow signal data to be inputed at NN
    
#     #read the first na values of Yr vector and returns a vector reshaped to data train
#     for i = 1:ns-1
#         H[1:na*ny,i] = reshape(Yr[:,(1:na).+(i-1)],(na*ny,1))
#         H[na*ny.+(1:nu*nb),i] = reshape(U[:,(1:nb).+(i-1)],(nb*nu,1))
#     end
#     H = [vec(H[:,i]) for i =1:ns-1] |> gpu

#     return H, Yt
# end

# function TrainNARX(model,Data)
#     opt_state = Flux.setup(Adam(), model.nn) 
#     Flux.train!(model.nn, Data, opt_state) do m, x, y
#         Flux.mse(m(x), y)
#     end
#     #return model.nn
# end

# function nnlsim(model, U, x0)
#     nu = model.nu
#     ny = model.ny
#     na = model.na;
#     nb = model.nb;
    
#     r1 = na-1;
#     r2 = nb-1
#     Ypred = zeros(ny,ns+r1) |> gpu;
#     Ugpu = zeros(nu,ns+r2) |> gpu;
#     CUDA.@allowscalar Ypred[:,r1] = x0;
    

#     for i = 1:ns-1
#         Yinp = [reshape(Ypred[:,(1:na).+(i-1)],(ny*na)); reshape(Ugpu[:,(1:nb).+(i-1)],(nu*nb))] |> gpu
#         CUDA.@allowscalar Ypred[:,i+na] = model.nn(Yinp);
#     end
#     return Ypred[:,(1:ns).+r1]
# end

# function OptimizeNARX(Y,U,na, nb,lim)
#     ny,ns = size(Y); #determines the number of outputs and the number of samples
#     nu,_ = size(U); #determines the number of inputs
#     #if the signal or the input come as rows being the samples, it changes the order for columns being samples
#     Y,ny,ns = CheckData(Y,ny,ns); 
#     U,nu,ns = CheckData(U,nu,ns);
#     H,Yt = DataNARX(Y,U,ny,nu,na,nb,ns)
#     Data = zip(H,Yt);
#     ho = @hyperopt for i=lim^2,
#         p = 1:lim

#         V = []
#         model, opt_state = NARX(ny,nu,na,nb,V);
#         Flux.train!(model, Data, opt_state) do m, x, y
#             Flux.mse(m(x), y)
#         end
#         cost = norm(ModelSimulation(model,ns,H)-Y)
#     end
#     p,q = ho.minimizer;

#     return TrainNARX(ny,nu,na,nb,p,q,Data)

# end

# function NNARX_LSTM(ny,nu,na,nb,p)
#     dataSize = ny*na+nu*nb
#     model = Chain(
#         LSTM(dataSize,dataSize), 
#         LSTM(dataSize,1), 
#         Dense(1,p,tanh),
#         Dense(p,2)
#         ) |> gpu;
#     return NNARXModel(model,na,nb,nu,ny)
# end

# function DataNorm(Y)
#     Ymax = maximum(abs.(Y),dims=2);
#     Ymin = minimum(abs.(Y),dims=2);
#     return (Y.-Ymin)./(Ymax-Ymin)

    
# end


# #-----------------------------KNN

# function KNNFilter(Yout,k)

#     kdtree = KDTree(Yout)
#     index_knn, _ = knn(kdtree,Yout,k,true)
#     Y_clean = zeros(ny,ns)
#     for i = 1:ns
#         Y_clean[i] = sum(Yout[index_knn[i]])/k
#     end
#     return Y_clean
# end

# function GRA_KNN(Y,U,nx,dt,p,k)
#     Yid = KNNFilter(Y,k)
#     return GRA(Yid,U,nx,dt,p);   
# end


# function ARXKD(Y,U, na, nb,dt)
#     ny,ns = size(Y); #determines the number of outputs and the number of samples
#     nu,_ = size(U); #determines the number of inputs
#     #if the signal or the input come as rows being the samples, it changes the order for columns being samples
#     Y,ny,ns = CheckData(Y,ny,ns); 
#     U,nu,ns = CheckData(U,nu,ns);

#     linhas = (ns-na)*ny; #since we want na coefficients to our denominator, we need to organize the samples.
#     colunas = na+nb*ny*nu; #it's related to the number of coefficients. 
#     H = zeros(linhas,colunas); #hankel matrix is reserved in memory

#     #for building our StateSpace system, we storage a matrix with the TF discrete type, to populate the array
#     G = Matrix{TransferFunction{Discrete{Float64}, ControlSystemsBase.SisoRational{Float64}}}(undef,ny,nu); 
    

#     #professor kurka formulation for MIMO-ARX with least square
#     Id = Matrix{Float64}(I,ny,ny);
#     H[:,1:nb*ny*nu] = kron(Hankel_PD2(U,nb,ns-na),Id);
#     H[:,(1:na).+nb*ny*nu] = -Hankel_PD(Y,na,linhas);

#     return KDTree(H')
# end

# # function ()
# #     x0 =vec([0 1 0])

# #     for i=1:ns
# #         index,_ = knn(kdtree,x0,k,true)
# #         Yesc = sum(Y[index])/k
# #         x0 = vec([0 x0[3] Yesc])

# #     end

    
# # end

