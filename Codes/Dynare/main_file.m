close all
clear all


alpha = 1/3;
beta = 0.99;
delta = 0.025;
rho = 0.95;
sigma = 1;
sigmae_vec =[0.01, 0.02, 0.03]; % Desviación estándar del shock


for i =1:3
    sigmae= sigmae_vec(i);
    save param_nc alpha beta delta rho sigma sigmae
    dynare Ejemplo_Dynare_2.mod

end