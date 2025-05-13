% Variables endógenas
var y I k a c w R r;

% Variables exógenas
varexo  e ;

% Parametros

parameters alpha beta delta rho sigma sigmae;

% Asignar valores a los parámetros

alpha = 1/3;
beta = 0.99;
delta = 0.025;
rho = 0.95;
sigma = 1;
sigmae = 0.01;

% Bloque del modelo

model;
exp(c)^(-sigma) = beta*exp(c(+1))^(-sigma)*(alpha*exp(a(+1))*exp(k)^(alpha-1) + (1-delta));
exp(k) = exp(a)*exp(k(-1))^(alpha) - exp(c) + (1-delta)*exp(k(-1));
a = rho*a(-1) + e;
exp(y) = exp(a)*exp(k(-1))^(alpha);
exp(I) = exp(y) - exp(c);
exp(c)^(-sigma) = beta*exp(c(+1))^(-sigma)*(1 + r);
exp(R) = alpha*exp(a)*exp(k(-1))^(alpha-1);
exp(w) = (1-alpha)*exp(a)*exp(k(-1))^(alpha);
end;

% Valores iniciales
initval;
k = log(30);
y = log(3);
c = log(2.5);
I = log(0.5);
a = 0;
r = (1/beta) - 1;
R = log((1/beta) - (1 - delta));
w = log(1);
end;

% Bloque Shocks
shocks;
var e = sigmae^2;
end;

% Estado Estacionario 
steady;

stoch_simul(order=1, irf=200, periods=20000);



