function [T_order, T] = dynamic_resid_tt(y, x, params, steady_state, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 6
    T = [T; NaN(6 - size(T, 1), 1)];
end
T(1) = exp(y(13))^(-params(5));
T(2) = params(2)*exp(y(21))^(-params(5));
T(3) = params(1)*exp(y(20))*exp(y(11))^(params(1)-1);
T(4) = exp(y(3))^params(1);
T(5) = exp(y(12))*T(4);
T(6) = params(1)*exp(y(12))*exp(y(3))^(params(1)-1);
end
