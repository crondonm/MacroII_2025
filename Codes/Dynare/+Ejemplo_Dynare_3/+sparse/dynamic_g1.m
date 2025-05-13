function [g1, T_order, T] = dynamic_g1(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T_order, T)
if nargin < 9
    T_order = -1;
    T = NaN(10, 1);
end
[T_order, T] = Ejemplo_Dynare_3.sparse.dynamic_g1_tt(y, x, params, steady_state, T_order, T);
g1_v = NaN(26, 1);
g1_v(1)=(-((1-params(3))*exp(y(3))+T(8)));
g1_v(2)=(-T(8));
g1_v(3)=(-(params(1)*exp(y(12))*exp(y(3))*getPowerDeriv(exp(y(3)),params(1)-1,1)));
g1_v(4)=(-(exp(y(12))*(1-params(1))*T(7)));
g1_v(5)=(-params(4));
g1_v(6)=exp(y(9));
g1_v(7)=(-exp(y(9)));
g1_v(8)=exp(y(10));
g1_v(9)=(-(T(2)*params(1)*exp(y(20))*exp(y(11))*getPowerDeriv(exp(y(11)),params(1)-1,1)));
g1_v(10)=exp(y(11));
g1_v(11)=(-T(5));
g1_v(12)=1;
g1_v(13)=(-T(5));
g1_v(14)=(-T(6));
g1_v(15)=(-(T(4)*exp(y(12))*(1-params(1))));
g1_v(16)=T(9);
g1_v(17)=exp(y(13));
g1_v(18)=exp(y(13));
g1_v(19)=T(9);
g1_v(20)=exp(y(14));
g1_v(21)=exp(y(15));
g1_v(22)=(-T(2));
g1_v(23)=(-(T(2)*T(3)));
g1_v(24)=(-((T(3)+1-params(3))*T(10)));
g1_v(25)=(-((1+y(16))*T(10)));
g1_v(26)=(-1);
if ~isoctave && matlab_ver_less_than('9.8')
    sparse_rowval = double(sparse_rowval);
    sparse_colval = double(sparse_colval);
end
g1 = sparse(sparse_rowval, sparse_colval, g1_v, 8, 25);
end
