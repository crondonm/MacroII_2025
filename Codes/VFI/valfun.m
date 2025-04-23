function val=valfun(k)

% Test 

% This program gets the value function for a neoclassical growth model with
% no uncertainty and CRRA utility

global v0 beta delta alpha kmat k0 sigma

g=interp1(kmat,v0,k,'linear'); % smooths out previous value funtion

% defining consumption

c=k0^alpha-k+(1-delta)*k0; % Note that k is the input

% For a given K, c can reach negative values, therefore we need to make
% sure this values are not going to be chosen during the maximization.

if c<=0
    val=-8888888888-800*abs(c);
else
    val=((1/(1-sigma))*((c^(1-sigma))-1))+beta*g;
    %val=log(c)+beta*g;
end

val=-val; % needs to be negative since the code minimizes this function

