mutable struct Controller
    M::Matrix
    L::Vector
    θ::Vector
    h::Vector
    na::Int64  
    nb::Int64
    η::Float64
    i::Int64
end

function ARX_oloop(y,u,Control)
    h = Control.h
    M = Control.m
    L = Control.L
    θ = Control.θ
    η = Control.η
    Control.i += 1
    i = Control.i

    #VERIFICAR A CADEIA DAS ESTRUTURAS h
    #VERIFICAR COERÊNCIA DO ALGORITMO NESTE SENTIDO
    


    h[1:na-1] = h[2:na];
    h[na] = y;
    h[na+nb] = u; #apenas para o caso onde nb =1

    M .+= h*h' #aqui é o inverso da formulação, pois vetores no julia são verticais
    L .+= h*y 
    ∇ = 2*(M*θ-L)/i
    θ = θ - η*∇

    
end

function EcoSystem(Y,U,na,nb,η)

    M = zeros(na+nb,na+nb);
    L = zeros(na+nb);
    θ = rand(na+nb);
    h = zeros(na+nb)
    i = 0
    η = 0.001
    Arduino = Controller(M,L,θ,h,na,nb,η,i)
    m = 0
    while true

        y = Y[m]
        u = 
        ARX_oloop(y,u,Arduino)
        m += 1
    end


    

end