var Y I K N A C W Rk r logC logY logI logK logN logW logA alp;

varexo eA;

parameters beta alpha delta rhoA SA ns ks theta ys is cs ws;

load param.rbc.moments;
set_param_value('alpha', alpha);
set_param_value('beta', beta);
set_param_value('delta', delta);
set_param_value('rhoA', rhoA);
set_param_value('SA', SA);
set_param_value('ns', ns);
set_param_value('ks', ks);
set_param_value('theta', theta);
set_param_value('ys', ys);
set_param_value('is', is);
set_param_value('cs', cs);
set_param_value('ws', ws);

model;

% (1) Euler equation bonds
(1/C) = beta*(1/(C(+1)))*(Rk(+1)+(1-delta));

% (2) Labor supply
theta/(1-N) = (1/C)*W; 

% (3) Labor demand
W = (1-alpha)*A*K(-1)^(alpha)*N^(-alpha); 

% (4) Capital demand
Rk = alpha*A*K(-1)^(alpha-1)*N^(1-alpha); 

% (5) Output
Y = A*K(-1)^(alpha)*N^(1-alpha);

%(6) Resource
Y = C+I; 

% (7) Law of motion capital
K = I+(1-delta)*K(-1); 

% (8) Euler equation bonds
(1/C) = beta*(1/C(+1))*(1+r);

% (9) TFP process
log(A) = rhoA*log(A(-1)) + SA*eA; 

% Define logs
logC = log(C); 
logY = log(Y); 
logI = log(I); 
logN = log(N);
logW = log(W); 
logA = log(A); 
logK = log(K); 

% average labor productivity
alp = logY-logN; 
end;

shocks;
var eA = 1; 
end;

initval;
Y = ys; 
C = cs; 
I = is; 
N = ns; 
K = ks; 
W = ws; 
r = (1/beta-1); 
Rk = (1/beta-1+delta); 
A = 1; 
logA = 0; 
logY = log(ys); 
logC = log(cs); 
logI = log(is); 
logK = log(ks); 
logW = log(ws); 
alp = log(ys)-log(ns); 

end;

steady;

stoch_simul(order=1, irf=20, nograph, ar=4, hp_filter=1600); 