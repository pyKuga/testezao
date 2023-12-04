mutable struct Controller
    M::Matrix
    L::Vector
    h::Vector
    θ::Vector
    ∇::Vector
    na::Int64  
    nb::Int64
    i::Int64
    η::Float64
    tol::Float64
end

mutable struct AdamOpt
    η::Float64
    α::Float64
    β::Float64
    t::Int64
    ϵ::Float64
    m::Vector
    v::Vector
end

function ControllerInit(η,na,nb,tol)
    z = zeros(na+nb);
    r = rand(na+nb)
    i = na
    return Controller(z*z',z,z,r,r,na,nb,i,η,tol)
end

function  MSEGrad(c,ns)
    return 2*(c.M*c.θ-c.L)/ns
end

function AdamInit(η,α,β,ϵ)
    t = 0
    z = zeros(na+nb)
    AdamOpt(η,α,β,t,ϵ,z,z)
end

function AdamRun(c)
    opt = AdamInit(c.η,0.9,0.999,1e-8)
    Δ = 100
    while norm(Δ)/norm(c.θ) > c.tol
        opt.t += 1;
        c.∇ = MSEGrad(c,ns);
        opt.m = opt.α*opt.m + (1-opt.α)*c.∇;
        opt.v = opt.β*opt.v + (1-opt.β)*c.∇.^2;
        hat_m = opt.m/(1-opt.α^opt.t);
        hat_v = opt.v/(1-opt.β^opt.t);
        Δ = c.η*hat_m./(sqrt.(hat_v).+opt.ϵ)
        c.θ = c.θ - Δ
        print("\r"*string(norm(Δ)/norm(c.θ)))
    end
end

function grad_Bulid(Arduino,ns,Y,U)
    h = Arduino.h;
    M = Arduino.M;
    L = Arduino.L;
    i = Arduino.i;

    h[1:na] = Y[1:na];

    while i < ns
        i += 1
        h[na+nb] = U[i]; #apenas para o caso onde nb =1
        M .+= h*h'; #aqui é o inverso da formulação, pois vetores no julia são verticais
        L .+= h*Y[i] ;
        h[1:na-1] = h[2:na];
        h[na] = Y[i];
    end
end
# L = zeros(na+nb)
# M = zeros(na+nb,na+nb)

function GDS(c)
    while norm(c.η*c.∇)/norm(c.θ) > c.tol
        c.∇ = MSEGrad(c,ns)
        c.θ = c.θ - c.η*c.∇
        print("\r"*string(norm(c.η*c.∇)/norm(c.θ)))
    end
end




