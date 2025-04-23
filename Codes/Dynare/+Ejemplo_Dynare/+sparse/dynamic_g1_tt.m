function [T_order, T] = dynamic_g1_tt(y, x, params, steady_state, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = Ejemplo_Dynare.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
T_order = 1;
if size(T, 1) < 10
    T = [T; NaN(10 - size(T, 1), 1)];
end
T(7) = exp(y(3))*getPowerDeriv(exp(y(3)),params(1),1);
T(8) = exp(y(12))*T(7);
T(9) = exp(y(13))*getPowerDeriv(exp(y(13)),(-params(5)),1);
T(10) = params(2)*exp(y(21))*getPowerDeriv(exp(y(21)),(-params(5)),1);
end
