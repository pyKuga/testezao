mutable struct Controller
    M::Matrix
    L::Vector
    θ::Vector
    h::Vector
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
    m_p::Vector
    m_c::Vector
    v_p::Vector
    v_c::Vector
end

function ControllerInit(η,na,nb)
    z = zeros(na+nb);
    r = rand(na+nb)
    i = na
    tol = 10e-8
    return Controller(
        z*z',
        z,
        r,
        z,
        r,
        na,
        nb,
        i,
        η,
        tol
        )
end

function  MSEGrad(c,ns)
    return 2*(c.M*c.θ-c.L)/ns
end

function AdamInit(na,nb,η,α,β)
    t = 0
    ϵ = 1e-8
    z = zeros(na+nb)
    AdamOpt(η,α,β,t,ϵ,z,z,z,z)
end

function AdamRun(c)
    opt = AdamInit(c.η,0.9,0.999,0,10e-8)
    while norm(c.Δ)/norm(c.θ) > c.tol
        t += 1;
        c.∇ = MSEGrad(Controller,ns);
        opt.m_c = opt.α*opt.m_p + (1-opt.α)*∇;
        opt.v_c = opt.β*opt.v_p + (1-opt.β)*∇.^2;
        hat_mc = opt.m_c/(1-α^t);
        hat_vc = opt.v_c/(1-β^t);
        Δ = c.η*hat_mc./(sqrt.(hat_vc).+ϵ)
        c.θ = c.θ - Δ
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
        c.∇ = MSEGrad(Controller,ns)
        c.θ = c.θ - c.η*c.∇
    end
end




