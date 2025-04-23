function [y, T, residual, g1] = static_2(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(2, 1);
  T(1)=exp(y(5));
  T(2)=exp(y(4));
  T(3)=exp(y(3));
  T(4)=T(3)^params(1);
  residual(1)=(T(3))-(T(2)*T(4)-T(1)+T(3)*(1-params(3))+x(2));
  T(5)=T(1)^(-params(5));
  T(6)=params(1)*T(2)*T(3)^(params(1)-1);
  T(7)=T(6)+1-params(3);
  residual(2)=(T(5))-(T(5)*params(2)*T(7));
  T(8)=T(1)*getPowerDeriv(T(1),(-params(5)),1);
if nargout > 3
    g1_v = NaN(4, 1);
g1_v(1)=T(1);
g1_v(2)=T(8)-T(7)*params(2)*T(8);
g1_v(3)=T(3)-(T(3)*(1-params(3))+T(2)*T(3)*getPowerDeriv(T(3),params(1),1));
g1_v(4)=(-(T(5)*params(2)*params(1)*T(2)*T(3)*getPowerDeriv(T(3),params(1)-1,1)));
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 2, 2);
end
end
