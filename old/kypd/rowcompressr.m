function [Ac,U,r] = rowcompressr(A,tol)
%
% [Ac,U,r] = rowcompressr(A,tol)
%
% Computes U such that U'*A = [0; Ac] with Ac having full row rank r
%

% Copyright (c) by Anders Hansson 1996

[m,n] = size(A);
[U,S,V] = svd(A);
r = rank(S,tol*eps);
U1 = U(:,1:r);
U2 = U(:,r+1:m);
Ac = U1'*A;
U = [U2 U1];