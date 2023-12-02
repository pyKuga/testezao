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

B = [0 ; 1]

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
sysARX = ARX_K(Y,U,na,nb,dt)

