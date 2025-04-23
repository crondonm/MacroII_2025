function [residual, T_order, T] = static_resid(y, x, params, T_order, T)
if nargin < 5
    T_order = -1;
    T = NaN(4, 1);
end
[T_order, T] = Ejemplo_Dynare.sparse.static_resid_tt(y, x, params, T_order, T);
residual = NaN(8, 1);
    residual(1) = (T(1)) - (T(1)*params(2)*(T(2)+1-params(3)));
    residual(2) = (exp(y(3))) - (T(4)-exp(y(5))+exp(y(3))*(1-params(3))+x(2));
    residual(3) = (y(4)) - (y(4)*params(4)+x(1));
    residual(4) = (exp(y(1))) - (T(4));
    residual(5) = (exp(y(2))) - (exp(y(1))-exp(y(5)));
    residual(6) = (T(1)) - (T(1)*params(2)*(1+y(8)));
    residual(7) = (exp(y(7))) - (T(2));
    residual(8) = (exp(y(6))) - (T(3)*exp(y(4))*(1-params(1)));
end
