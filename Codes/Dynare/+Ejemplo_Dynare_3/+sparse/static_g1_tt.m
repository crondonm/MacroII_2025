function [T_order, T] = static_g1_tt(y, x, params, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = Ejemplo_Dynare_3.sparse.static_resid_tt(y, x, params, T_order, T);
T_order = 1;
if size(T, 1) < 8
    T = [T; NaN(8 - size(T, 1), 1)];
end
T(5) = params(1)*exp(y(4))*exp(y(3))*getPowerDeriv(exp(y(3)),params(1)-1,1);
T(6) = exp(y(3))*getPowerDeriv(exp(y(3)),params(1),1);
T(7) = exp(y(4))*T(6);
T(8) = exp(y(5))*getPowerDeriv(exp(y(5)),(-params(5)),1);
end
