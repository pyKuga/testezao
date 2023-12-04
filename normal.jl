using LinearAlgebra
using ControlSystems
using Plots
using Random
using Statistics

include("sysid.jl")
include("KugaPack.jl")
include("ARX_ON.jl")

m = 1
c = 2
k = 50

A = [
    0 1;
    -k/m -c/m
]

B = [0 ; 1/m]

C = [1 0]

D = 0

sys = ss(A,B,C,D)

#numerical conditions are set up
dt = 2e-3;
maxT = 10;
ns = round(Int,maxT/dt)+1;
t = 0:dt:maxT;

#problem conditions are set up
x0 = [0; 0];

#identification parameters
U = 4*PulseGen(ns); 
Y,_,_ = lsim(sys,U,t,x0);

na= 2
nb = 1
η = 0.001
tol = 1e-8
sysARX = ARX_K(Y,U,na,nb,dt);

Arduino = ControllerInit(η,na,nb);

grad_Bulid(Arduino,ns,Y,U)

M =  [40000.0    -787.791   -788.741; -787.791    20.5128    20.5123; -788.741    20.5123    20.5128]
L =  [789.6883998477451 ,-20.510924087846533,-20.512348335936757]

Arduino.M = M
Arduino.L = L

GDS(c)
AdamRun(c)



# Arduino.M\Arduino.L

# Identified = tf(Arduino.θ[na+nb],[1; -Arduino.θ[1:na]],dt)
# dampreport(Identified)


# Ys,_,_ = lsim(Identified,U,t,x0);
# plot(t,Ys')


