function [residual, T_order, T] = dynamic_resid(y, x, params, steady_state, T_order, T)
if nargin < 6
    T_order = -1;
    T = NaN(6, 1);
end
[T_order, T] = Ejemplo_Dynare_3.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
residual = NaN(8, 1);
    residual(1) = (T(1)) - (T(2)*(T(3)+1-params(3)));
    residual(2) = (exp(y(11))) - (T(5)-exp(y(13))+(1-params(3))*exp(y(3)));
    residual(3) = (y(12)) - (params(4)*y(4)+x(1));
    residual(4) = (exp(y(9))) - (T(5));
    residual(5) = (exp(y(10))) - (exp(y(9))-exp(y(13)));
    residual(6) = (T(1)) - (T(2)*(1+y(16)));
    residual(7) = (exp(y(15))) - (T(6));
    residual(8) = (exp(y(14))) - (T(4)*exp(y(12))*(1-params(1)));
end
