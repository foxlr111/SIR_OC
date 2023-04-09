%% Simulation of the SIR Model

clear all
close all
clc


%% Initialize 

% time
dt = 0.05;          
tvals = 0:dt:150;  

% parameters
beta = 0.00009;
gamma = 0.05;

pars = [beta; gamma];

% initial conditions
N = 10000;
S0 = N-1;
I0 = 1;
R0 = 0;
Cases0 = 0;

inits = [S0; I0; R0; Cases0];


%% Solve Model

[time, ODEmodel] = ode45(@(t,x)SIR_Model(t,x,pars), tvals, inits);
S = ODEmodel(:,1);
I = ODEmodel(:,2);
R = ODEmodel(:,3);
Cases = ODEmodel(:,4);


%% Total Cases

Cases(end)

% peak
max(I)
find(I==max(I))*dt


%% Plot

hold on
P1 = plot(time,S);
set(P1, "linewidth",2,"color","#2623D6") 
P2 = plot(time,I);
set(P2, "linewidth",2,"color","#D62323") 
P3 = plot(time,R);
set(P3,"linewidth",2,"color","#468220") 
set(gca,"fontsize",12)
title("SIR model","fontsize",20)
xlabel("days","fontsize",16)
ylabel("individuals","fontsize",16)
LEG = legend("S","I","R");
set(LEG,"fontsize",14)
hold off


%% Model

function dx = SIR_Model(t,x,pars)

% parameters
beta = pars(1);
gamma = pars(2);

% states
S = x(1);
I = x(2);
R = x(3);
Cases = x(4);

% equations
dS = -beta*I*S;
dI = beta*I*S - gamma*I;
dR = gamma*I;
dCases = beta*I*S;
dx = [dS; dI; dR; dCases];
end