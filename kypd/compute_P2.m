function P2=compute_P2(P1,A,B,c,M,solver)
  
% Computes the part of P corresponding to the uncontrollable part
% of A. A is partitiond as A=[A1  A2] and B as B=[B1].
%                            [ 0  A2]            [ 0]
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n,m]=size(B);

A1=A(1:c,1:c);
A12=A(1:c,c+1:n);
A2=A(c+1:n,c+1:n);
B1=B(1:c,:);

Q1=M(1:c,1:c);
Q12=M(1:c,c+1:n);
Q2=M(c+1:n,c+1:n);
S1=M(1:c,n+1:n+m);
S2=M(c+1:n,n+1:n+m);
R=M(n+1:n+m,n+1:n+m);

temp=P1*B1+S1;
X=[A1'*P1+P1*A1+Q1 temp;temp' R];
X=0.5*(X+X');
Y=[P1*A12+Q12;S2'];
[ma,na] = size(A2);
[mb,nb] = size(A2');
if solver=='schur'
  [ua,ta] = schur(A2','complex');
  j = ma:-1:1;
  ub = ua(:,j);
  tb = ta(j,j)';
  P2=lyapunov(ua,ta,ub,tb,Q2-Y'*inv(X)*Y);
end;
if solver=='diago' | solver=='block'
  [ua,ta] = eig(A2'); 
  uainv = inv(ua);
  P2=diaglyapunov(ua,uainv,ta,Q2-Y'*inv(X)*Y);
end;