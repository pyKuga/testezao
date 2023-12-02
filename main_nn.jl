#School of Mechanical Engineering - State University of Campinas
#Paulo Yoshio Kuga
#First Release: October 2023

#Present version: 20231001

#this is the main program to run neural network identification 

using LinearAlgebra
using ControlSystems
using Plots
using ControlSystemIdentification
using Random
using Statistics
using DelimitedFiles
using CUDA
using Flux
using Hyperopt


include("KugaPack.jl")
include("sysid.jl")
include("blockpack.jl")
include("analysis.jl")


#numerical conditions are set up
dt = 1e-3;
maxT = 6;
ns = round(Int,maxT/dt)+1;
t = 0:dt:maxT;

#problem conditions are set up
x0 = [0; 0; 0; 0;];

w = 5
U = 4*ones(1,ns);
#ones(1,ns)'
#chirp(w,t)';
#0.01*PulseGen(ns); 

#using the given parameters, the system (discrete space state) is created 
Param = [0.77 0.59 2.1 1.2 9 200 200 200];
sysP,sysD = F2DOF(dt,Param);

Y,_,_= lsim(sysD,U,t,x0=x0); #time response
magD,phaseD,w = bode(sysD); #frequency response

#plot(t,Y')

#Since F2DOF is a two degree freedom problem, we can state that the number of outputs is 2
#This implies the noise beign as (2,n)
#In this analysis, we are using MersenneTwister as the random number generator. 

seed = 98832+5556594 #seed for mersenne twister
noise = noise_gen(seed,2,ns);

fineza = 1;

NAmp = [i*1e-2 for i in 0:fineza:100];

#identification parameters
nx = 4

na = 100;
nb = 1;
ny = 2;
nu = 1;
p = 25

x0 = [0; 0]
H, Yt = DataNARX(Y,U, ny, nu, na, nb,ns);
Data = zip(H,Yt);

model = NNARX_LSTM(ny,nu,na,nb,p)

opt_state = Flux.setup(Adam(), model.nn); 
Flux.train!(model.nn, Data, opt_state) do m, x, y
        Flux.mse(m(x), y);
end

Ypred = nnlsim(model, U, x0)

plot(t,[Y' CUDA.@allowscalar Ypred'])



Ygpu = Y |> gpu

cost = zeros(100)

for i=100:100:p
        model = NNARX_LSTM(ny,nu,na,nb,i)
        opt_state = Flux.setup(Adam(), model.nn); 
        Flux.train!(model.nn, Data, opt_state) do m, x, y
                Flux.mse(m(x), y);
        end
        cost[Int(i/100)] = norm(nnlsim(model, U, x0)-Ygpu);
        print(string(i) * "|")
end


plot(100:100:p,cost)

CSV.write("100to10000.csv", cost)
writedlm( "100to10000.csv",  cost, ',')



for pair in Data
        # Unpack this element (for supervised training):
        h, yr = pair
      
        # Calculate the gradient of the objective
        # with respect to the parameters within the model:
        grads = Flux.gradient(model) do m
            result = m(h)
            sGRA = GRA(Yout,U,nx,dt,p);
            loss(h, yr)
        end
      
        # Update the parameters so as to reduce the objective,
        # according the chosen optimisation rule:
        Flux.update!(opt_state, model, grads[1])
      end