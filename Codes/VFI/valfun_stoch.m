function val=valfun_stoch(k)

%  Value Function Iteration: Stochastic growth model with CRRA utility

global v0 beta delta alpha kmat prob a0 j sigma k0

g=interp1(kmat,v0,k,'linear'); % smooths out previous value function

c=a0*k0^alpha-k+(1-delta)*k0; % consumption

if c<=0
    val=-8888888888888888-800*abs(c); % avoid negative consumption values
else
    val=(1/(1-sigma))*((c^(1-sigma))-1)+beta*(g*prob(j,:)') ;% Value function with the expected value included
    %val=log(c)+beta*(g*prob(j,:)') ;% Value function with the expected value included
end

val=-val; % negative since we are minimizing

