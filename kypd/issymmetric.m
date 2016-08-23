function f=issymmetric(M)

[n,m]=size(M);

f=0;

if m~=n
  error('The matrix is not square')
elseif norm(M-M',1)<n*eps
  f=1;
end;
