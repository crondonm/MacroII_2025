function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = Ejemplo_Dynare_3.static_resid_tt(T, y, x, params);
end
residual = zeros(8, 1);
    residual(1) = (T(1)) - (T(1)*params(2)*(T(2)+1-params(3)));
    residual(2) = (exp(y(3))) - (T(4)-exp(y(5))+exp(y(3))*(1-params(3)));
    residual(3) = (y(4)) - (y(4)*params(4)+x(1));
    residual(4) = (exp(y(1))) - (T(4));
    residual(5) = (exp(y(2))) - (exp(y(1))-exp(y(5)));
    residual(6) = (T(1)) - (T(1)*params(2)*(1+y(8)));
    residual(7) = (exp(y(7))) - (T(2));
    residual(8) = (exp(y(6))) - (T(3)*exp(y(4))*(1-params(1)));

end
