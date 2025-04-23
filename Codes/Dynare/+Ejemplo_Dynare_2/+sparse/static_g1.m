function [g1, T_order, T] = static_g1(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T_order, T)
if nargin < 8
    T_order = -1;
    T = NaN(8, 1);
end
[T_order, T] = Ejemplo_Dynare_2.sparse.static_g1_tt(y, x, params, T_order, T);
g1_v = NaN(21, 1);
g1_v(1)=exp(y(1));
g1_v(2)=(-exp(y(1)));
g1_v(3)=exp(y(2));
g1_v(4)=(-(T(1)*params(2)*T(5)));
g1_v(5)=exp(y(3))-(exp(y(3))*(1-params(3))+T(7));
g1_v(6)=(-T(7));
g1_v(7)=(-T(5));
g1_v(8)=(-(exp(y(4))*(1-params(1))*T(6)));
g1_v(9)=(-(T(1)*params(2)*T(2)));
g1_v(10)=(-T(4));
g1_v(11)=1-params(4);
g1_v(12)=(-T(4));
g1_v(13)=(-T(2));
g1_v(14)=(-(T(3)*exp(y(4))*(1-params(1))));
g1_v(15)=T(8)-(T(2)+1-params(3))*params(2)*T(8);
g1_v(16)=exp(y(5));
g1_v(17)=exp(y(5));
g1_v(18)=T(8)-(1+y(8))*params(2)*T(8);
g1_v(19)=exp(y(6));
g1_v(20)=exp(y(7));
g1_v(21)=(-(T(1)*params(2)));
if ~isoctave && matlab_ver_less_than('9.8')
    sparse_rowval = double(sparse_rowval);
    sparse_colval = double(sparse_colval);
end
g1 = sparse(sparse_rowval, sparse_colval, g1_v, 8, 8);
end
