function [check,matrix_info,T,c]=check_A(matrix_info,i)

% [check,matrix_info,T,c]=check_A(matrix_info,i)
% 
% Checks if the system is controllable or stabilizable. If check=0
% the system is not stabilizable, if check=1 the system is stabilizable
% and if check=2 the system is controllable. If the system is only
% stabilizable the system is transformed.
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
K=matrix_info.K;
n=matrix_info.n(i);
nm=matrix_info.nm(i);
m=nm-n;
A=matrix_info.A{i};
B=matrix_info.B{i};

[A,B,T,c] = vandooren1(A,B,10);
A=fliplr(flipud(A));
B=flipud(B);
T=fliplr(T);
if c<n&min(real(eig(A(c+1:n,c+1:n))))>0 
    check=0;
elseif c<n
    check=1;
    matrix_info.A{i}=A;
    matrix_info.B{i}=B;
    Tbar=blkdiag(T,eye(m));
    matrix_info.M0{i}=Tbar'*matrix_info.M0{i}*Tbar;
    for j=1:K
        matrix_info.M{i,j}=Tbar'*matrix_info.M{i,j}*Tbar;
    end;
    matrix_info.C{i} = T*matrix_info.C{i}*T';
else
    check=2;
    T=eye(n);
end;


function [Ac,Bc,T,c] = vandooren1(A,B,tol)
%
% [Ac,Bc,T,c] = vandooren1(A,B,tol)
%
% Computes a state transformation T such that 
%
% A*T = T*Ac; B = T*Bc, where 
%
% Ac = [Ac1  0  ]; Bc = [0  ] 
%      [Ac21 Ac2]       [Bc2]
%
% with (Ac2,Bc2) controllable and where the dimension of Ac2 is c
%
% The routine implemented is described in van Dooren, Paul M.: The Generalized
% Eigenstructure Problem in Linear System Theory, IEEE Transaction on Automatic
% Control, Vol. AC-26, No. 1, February 1981
%

% Copyright (c) by Anders Hansson 1996

[n,m] = size(B);
[n,dummy] = size(A);
c = 0;
tau = n;
rho = m;
T = eye(n);
P = [A B];
ready = 0;
while ~ready,
   [barB,U,newrho] = rowcompressr(P(:,tau+1:tau+rho),tol);
   newtau = tau-newrho;
   if newrho == 0,
      c = n-tau;
      ready = 1;
   elseif newtau == 0,
      c = n;
      ready =1;
   else
      P = U'*P(:,1:tau)*U;
      P = P(1:newtau,:);
      T = T*[U zeros(tau,c); zeros(c,tau) eye(c)];
      tau = newtau;
      rho = newrho;
      c = c+rho;
   end
end
Ac = T'*A*T;
Bc = T'*B;