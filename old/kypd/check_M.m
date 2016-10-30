function [check,matrix_info,Pbar,V]=check_M(matrix_info,solver,lowrank,tol)
 
% [check,matrix_info,Pbar,V]=check_M(matrix_info,tol)
%
% Eliminates the (1,1)-block of the M-matrices and checks if they
% are linearly independent with tolerance tol. If not the number of
% M-matrices are reduced if possible. The values of check can be 0
% which means that the matrices are linearly dependent but the
% system can not be reduced, 1 which means that the matrices are
% linearly dependent and the system is reduced or 2 which means
% that the matrices are linearly independent.
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
N=matrix_info.N;
K=matrix_info.K;
nm=matrix_info.nm;
n=matrix_info.n;
m=nm-n;

if solver == 'schur' & ~lowrank
    matrix_info.schur.ua = cell(N,1);
    matrix_info.schur.ta = cell(N,1);
    matrix_info.schur.ub = cell(N,1);
    matrix_info.schur.tb = cell(N,1);
elseif solver == 'diago' & ~lowrank
    matrix_info.diago.ua = cell(N,1);
    matrix_info.diago.ta = cell(N,1);
    matrix_info.diago.uainv = cell(N,1);
elseif lowrank
    matrix_info.block_info = zeros(N,1);
    matrix_info.diago.ta = cell(N,1);
end   
      
matrix_info.T_diag = cell(N,1);
matrix_info.T_diaginv = cell(N,1);
matrix_info.D_diag = cell(N,1);

check=2;

% Eliminate the (1,1)-block of the M-matrices.

Pbar=[];
vec_M=[];
D_diag = cell(N,1);
c = matrix_info.c; %JANNE

for i=1:N
  temp=[];
  M0=matrix_info.M0{i};
  
  if n(i)~=0 
    A=matrix_info.A{i};
    B=matrix_info.B{i};
    [ma,na] = size(A);
    [mb,nb] = size(A);
    
    if solver=='schur'
      [ua,ta] = schur(A','complex');
      j = ma:-1:1;
      ub = ua(:,j);
      tb = ta(j,j)';
      Pbar{i,1}=lyapunov(ua,ta,ub,tb,-M0(1:n(i),1:n(i)));
      
      if lowrank
          [T,D] = eig(A');
      end
    end;
    
    if solver=='diago'
      [ua,ta] = eig(A');
      uainv = inv(ua);
      Pbar{i,1}=diaglyapunov(ua,uainv,ta,-M0(1:n(i),1:n(i)));
      
      if lowrank
          T = ua;
          D = ta;
      end
    end;
    
    if lowrank
        [matrix_info.T_diag{i},matrix_info.T_diaginv{i}, ...
            matrix_info.D_diag{i},matrix_info.block_info(i)] = create_blockdiag_sort(T,D);
    end
    
    M0(1:n(i),1:n(i))=zeros(n(i));
    M0(1:n(i),n(i)+1:nm(i))=M0(1:n(i),n(i)+1:nm(i))-Pbar{i,1}*B;
    M0(n(i)+1:nm(i),1:n(i))=M0(1:n(i),n(i)+1:nm(i))';
    matrix_info.M0{i}=M0;
  end;
  
    for j=1:K
        M=matrix_info.M{i,j};
        if n(i)~=0
            if solver=='schur'
                Pbar{i,j+1}=lyapunov(ua,ta,ub,tb,-M(1:n(i),1:n(i)));
            end;
            if solver=='diago'
            Pbar{i,j+1}=diaglyapunov(ua,uainv,ta,-M(1:n(i),1:n(i)));
            end;
            M(1:n(i),1:n(i))=zeros(n(i));
            M(1:n(i),n(i)+1:nm(i))=M(1:n(i),n(i)+1:nm(i))-Pbar{i,j+1}*B;
            M(n(i)+1:nm(i),1:nm(i))=M(1:nm(i),n(i)+1:nm(i))';
            matrix_info.M{i,j}=M;
        end;
        temp=[temp [vec(M(n(i)+1:nm(i),1:n(i)));...
		  hvec(M(n(i)+1:nm(i),n(i)+1:nm(i)))]];
      try
          C = matrix_info.C{i}; %JANNE
          c(j) = c(j) - trace(C*Pbar{i,j+1}); %JANNE
      catch
          %No C! c(j) needs not to be changed.
      end
    end; 
    
    vec_M=[vec_M;temp];
end;
matrix_info.c = c;
[n1,m1]=size(vec_M);

% Make a singular value decomposition.

[u,s,v]=svd(full(vec_M));
V=v';
if nargin<4
    tol = max(size(vec_M)')*max(max(s))*eps;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if the matrices are linearly independent.

r = sum(sum(s)>tol);
if r<m1
    disp('The M matrices are linearly dependent')
    
    % Check if the number of matrices can be reduced.
    
    if any(abs(V(r+1:m1,:)*matrix_info.c)>tol)
        check=0;
        error('The system cannot be reduced')
    else
      
        % Reduce the number of matrices.
	
        check=1;
        matrix_info.c=V(1:r,:)*matrix_info.c;
	matrix_info.K=r;
	E=u*s;
	E=E(:,1:r);
	% Create the new M matrices.
        for i=1:N
	    for j=1:r
	      matrix_info.M{i,j}=zeros(nm(i));
	      matrix_info.M{i,j}(n(i)+1:nm(i),1:n(i))=reshape(E(1: ...
	                                        n(i)*m(i),j),m(i),n(i));
	      matrix_info.M{i,j}=matrix_info.M{i,j}+matrix_info.M{i,j}';
	      matrix_info.M{i,j}(n(i)+1:nm(i),n(i)+1:nm(i))=...
		  inv_hvec(E(n(i)*m(i)+1:n(i)*m(i)+m(i)*(m(i)+1)/2,j));
	    end;
	    E=E(n(i)*m(i)+m(i)*(m(i)+1)/2+1:end,:);
	end;
    end;
end;

V=V(1:r,:);




function v=hvec(M)

[n,m]=size(M);

if n~=m
  error('The matrix is not square')
end;

if issymmetric(M)
  x=find(triu(ones(m)));
  v=M(x);
else
  error('The matrix is not symmetric')
end;




function M=inv_hvec(v)
  
m=(-1+sqrt(1+8*length(v)))/2;
if abs(round(m)-m)>1e-5
  error('The size of the vector is not compatible with a symmetric matrix')
else
  x=find(triu(ones(m)));
  M=zeros(m);
  M(x)=v;
  M=M+M'-diag(diag(M));
end;
