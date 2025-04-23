% Value Function Iteration
% Deterministic Growth Model

clear all
close all

% Initial Parameters

global v0 beta delta alpha kmat k0 sigma

sigma=2;
beta=0.95;
delta=0.1;
alpha=(1/3);

n=300; % # of observations in the grid for capital

% Solving by hand the k and c of the steady state

kstar=((alpha/((1/beta)-(1-delta))))^(1/(1-alpha));
cstar=(kstar^alpha)-delta*kstar;
ystar=kstar^alpha;

%% Creating a grid for capital 

kmin=0.25*kstar; % minimum value of capital
kmax=1.75*kstar; % maximum value of capital
stepsize=(kmax-kmin)/(n-1);
kmat=kmin:stepsize:kmax; % grid for capital

% kmat = linspace(kmin,kmax,n); % Matlab code

kmat=kmat'; % Column vector for capital n x 1

%% Value function iteration

% initial guess for value function and policy function

v0=zeros(n,1);
polfun=zeros(n,2);

% Parameters for iteration

tol=0.01; % tolerance
maxits= 1000; % Max number of iterations
dif=10000;
its=0;

while dif>tol & its<maxits
    for i=1:n
        k0=kmat(i,1);
        k1 = fminbnd(@valfun,kmin,kmax); % Value of Capital that maximizes consumption
        v1(i,1)=-valfun(k1); % Value of consumption given the capital we found in the last row
        k11(i,1)=k1; 
    end
 dif=norm(v1-v0);
 v0=v1;
 its=its+1;
end

 
for i = 1:n
   
    polfun(i,1) = kmat(i,1)^(alpha) - k11(i,1) + (1-delta)*kmat(i,1); % Consumption
    polfun(i,2) = k11(i,1); % K_{t+1}
    
end

figure
subplot(3,1,1)
plot(kmat,v0)
title('Value Function','Interpreter','latex')
subplot(3,1,2)
plot(kmat,polfun(:,1))
title('Consumption','Interpreter','latex')
subplot(3,1,3)
plot(kmat,k11)
title('Policy Function','Interpreter','latex')
suptitle('VFI: Deterministic Growth Model');

