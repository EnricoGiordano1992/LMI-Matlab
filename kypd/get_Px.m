function [u,P,x]=get_Px(X,F,G,matrix_info,n_basis)
 
% u is the objective function with the x recovered. 

N=matrix_info.N;
n=matrix_info.n;

b=[];
t=0;

for i=1:N
  for k=1:n_basis(i)
    b=[b;ip(X{i}-matrix_info.M0{i},F{t+k})];
  end;
  t=t+k;
end;

x=G\b;

tmp = 0;
for i=1:N
  if n(i)~0;
    sum=-X{i}+matrix_info.M0{i}+arraysum1(matrix_info,i,x);
    P{i}=lyap(matrix_info.A{i}',sum(1:n(i),1:n(i)));
    tmp = tmp+trace(matrix_info.C{i}*P{i});
  else
    P{i}=[];
  end;
end;


u=matrix_info.c'*x + tmp;
clear tmp