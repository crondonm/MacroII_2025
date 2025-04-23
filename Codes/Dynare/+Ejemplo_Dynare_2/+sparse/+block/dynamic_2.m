function [y, T, residual, g1] = dynamic_2(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(2, 1);
  T(1)=exp(y(13));
  T(2)=T(1)^(-params(5));
  T(3)=exp(y(21));
  T(4)=params(2)*T(3)^(-params(5));
  T(5)=params(1)*exp(y(20));
  T(6)=exp(y(11));
  T(7)=T(5)*T(6)^(params(1)-1)+1-params(3);
  residual(1)=(T(2))-(T(4)*T(7));
  T(8)=exp(y(12));
  T(9)=exp(y(3));
  T(10)=T(9)^params(1);
  residual(2)=(T(6))-(T(8)*T(10)-T(1)+(1-params(3))*T(9));
if nargout > 3
    g1_v = NaN(6, 1);
g1_v(1)=(-((1-params(3))*T(9)+T(8)*T(9)*getPowerDeriv(T(9),params(1),1)));
g1_v(2)=(-(T(4)*T(5)*T(6)*getPowerDeriv(T(6),params(1)-1,1)));
g1_v(3)=T(6);
g1_v(4)=T(1)*getPowerDeriv(T(1),(-params(5)),1);
g1_v(5)=T(1);
g1_v(6)=(-(T(7)*params(2)*T(3)*getPowerDeriv(T(3),(-params(5)),1)));
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 2, 6);
end
end
