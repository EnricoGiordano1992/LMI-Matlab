function [u,P,x,Z,soltime,errorflag]=kypd(matrix_info,options) 
 
if nargin<2
    options = sdpsettings;
end

for i=1:matrix_info.N
  matrix_info.M0{i}=matrix_info.M0{i}-options.kypd.tol*...
                    eye(size(matrix_info.M0{i})); 
end;

[F,n_basis,ntot_basis,F0]=basis_matrices(matrix_info,...
                                      options.kypd.lyapunovsolver,options.kypd.lowrank);
                                  
G=computeG(matrix_info,F,n_basis);
G0=computeG0(matrix_info,F0,ones(matrix_info.N,1));

t=0;
for i=1:matrix_info.N
  for j=1:n_basis(i)
    d(t+j,1)=ip(F{t+j},matrix_info.M0{i});
  end;
  t=t+n_basis(i);
end;

[X,Z,soltime,errorflag]=dual(d,F,G,n_basis,ntot_basis,matrix_info,options,F0,G0);

if ~errorflag
    [u,P,x]=get_Px(X,F,G,matrix_info,n_basis);
else 
    u = [];
    P = [];
    x = [];
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G=computeG(matrix_info,F,n_basis);

N=matrix_info.N;
K=matrix_info.K;

G=[];
Gtemp=[];
t=0;
for i=1:N
  for k=1:n_basis(i)
    for l=1:K
      Gtemp(k,l)=ip(F{t+k},matrix_info.M{i,l});
    end;
  end;
  t=t+k;
  G=[G;Gtemp];
  Gtemp=[];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G=computeG0(matrix_info,F0,n_basis);

     N=matrix_info.N;
     K=matrix_info.K;
 
 G=zeros(1,K);
 Gtemp=[];
 for i=1:N
     for l=1:K
       Gtemp(1,l)=ip(F0{i,1},matrix_info.M{i,l});
     end;
   G= G + Gtemp;
   Gtemp=[];
 end;
