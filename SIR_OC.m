%% Simulation of the SIR Model with Optimal Vaccination Rate

clear all
close all
clc


%% Initialize

% time          
dt = 0.05; 
T = 150;
tvals = 0:dt:T;       
tlen = length(tvals); 

% parameters
beta = 0.00009;
gamma = 0.05;
C1 = 1;
C2 = 0.1;
epsilon = 100;
vmax = 0.05;

pars = [beta; gamma; C1; C2; epsilon; vmax];

% initial conditions
N  = 10000;
S0 = N-1;
I0 = 1;
R0 = 0;
Cases0 = 0;

inits = [S0; I0; R0; Cases0];

% empty vectors
x = zeros(tlen,3); 
lambda = zeros(tlen,3);
v = zeros(tlen,1);


%% Forwards-Backwards Sweep

delta = 0.01;
test = -1;
count = 0;

while(test < 0 && count < 100)

    % set previous control, state, adjoint
    xold = x;
    lambdaold = lambda;
    vold = v; 
    
    % solve states (x)
    solx = ode45(@(t,x) states(t,x,pars,v,tvals), tvals, inits);
    x = deval(solx,tvals)';
    
    S = x(:,1);
    I = x(:,2);
    R = x(:,3);
    Cases = x(:,4);

    % solve adjoints (lambda)
    sollambda = ode45(@(t,lambda) adjoints(t,lambda,x,pars,v,tvals), [T 0], zeros(1,3));
    lambda = deval(sollambda,tvals)';
    
    lambdaS = lambda(:,1);
    lambdaI = lambda(:,2);
    lambdaR = lambda(:,3);
    
    % calculate control characteristic and update control  
    vchar = (-C2*S + lambdaS.*S - lambdaR.*S)/(2*epsilon);
    vstar = min(vmax, max(0, vchar));
    v = (vold + vstar)*0.5;
    
    % test convergence
    test = delta*sum(abs(v)) - sum(abs(vold - v));
    
    % update count
    count = count + 1
end


%% Total Cases

Cases(end)

% peak
max(I)
find(I==max(I))*dt


%% Total Cost (J)

Jcases = sum(C1*beta*I.*S).*dt
Jvax = sum(C2*v.*S + epsilon*v.^2).*dt
Jtot = Jcases + Jvax


%% Plots

subplot(2,1,1);
hold on
P1 = plot(tvals,S);
set(P1, "linewidth",2,"color","#2623D6") 
P2 = plot(tvals,I);
set(P2, "linewidth",2,"color","#D62323") 
P3 = plot(tvals,R);
set(P3,"linewidth",2,"color","#468220") 
set(gca,"fontsize",12)
title("(a) SIR model with vax","fontsize",20)
xlabel("days","fontsize",16)
ylabel("individuals","fontsize",16)
LEG = legend("S","I","R");
set(LEG,"fontsize",14)
hold off

subplot(2,1,2);
plot(tvals,v,"linewidth",2,"color","#9823D6")
set(gca,"fontsize",12)
title("(b) Optimal vax rate","fontsize",20)
xlabel("days","fontsize",16)
ylabel("prop. of ind. per day","fontsize",16)
LEG = legend("v");
set(LEG,"fontsize",14)


%% States

function dx = states(t,x,pars,v,tvals)

% parameters
beta = pars(1);
gamma = pars(2);
C1 = pars(3);
C2 = pars(4);
epsilon = pars(5);
vmax = pars(6);

% initialize
v = pchip(tvals,v,t);
dx = zeros(3,1);

% states
S = x(1);
I = x(2);
R = x(3);
Cases = x(4);

% state equations
dS = -beta*I*S - v*S;
dI = beta*I*S - gamma*I;
dR = gamma*I + v*S;
dCases = beta*I*S;
dx = [dS; dI; dR; dCases]; 
end


%% Adjoints

function dlambda = adjoints(t,lambda,x,pars,v,tvals)

% parameters
beta = pars(1);
gamma = pars(2);
C1 = pars(3);
C2 = pars(4);
epsilon = pars(5);
vmax = pars(6);

% initialize
dlambda = zeros(3,1);
x = interp1(tvals,x,t);
v = pchip(tvals,v,t);

% states
S = x(1);
I = x(2);
R = x(3);

% adjoints
lambdaS = lambda(1);
lambdaI = lambda(2);
lambdaR = lambda(3);

% adjoint equations
dlambdaS = -(C1*beta*I + C2*v - lambdaS*beta*I - lambdaS*v + lambdaI*beta*I + lambdaR*v);
dlambdaI = -(C1*beta*S - lambdaS*beta*S + lambdaI*beta*S - lambdaI*gamma + lambdaR*gamma);
dlambdaR = -(0);
dlambda  = [dlambdaS; dlambdaI; dlambdaR]; 
end
