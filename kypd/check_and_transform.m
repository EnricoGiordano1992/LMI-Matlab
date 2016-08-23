function [matrix_info,temp,T,Pbar,V,c,L,checkM,checkstab]=...
                   check_and_transform(matrix_info,transform,rho,solver,lowrank)
  
% [matrix_info,temp,T,Pbar,V,c,checkM,checkstab]=...
%                check_and_transform(matrix_info,transform,rho,solver)
%
% Checks if the M matrices are linearly independent, if the system
% is stabilizable or controllable and performs state feedback if
% necessary. If transrorm=1 state feedback is applied even if it is
% not necessary.
  
N=matrix_info.N;
K=matrix_info.K;
n=matrix_info.n;
nm=matrix_info.nm;
m=nm-n;

T{1}=[];
L{1}=[];
U{1}=[];
S{1}=[];

temp=matrix_info;

checkstab=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if the system is stabilizable or controllable.

for i =1:N
    if n(i)~=0
        [checkstab(i),matrix_info,T{i},c(i)]=check_A(matrix_info,i);
        if checkstab(i)==0
            error('The system is not stabilizable')
        elseif checkstab(i)==1 
            disp('The system is not controllable')
	    temp=matrix_info;
            matrix_info.A{i}=matrix_info.A{i}(1:c(i),1:c(i));
            matrix_info.B{i}=matrix_info.B{i}(1:c(i),:);
	    
	    matrix_info.M0{i}=matrix_info.M0{i}([1:c(i) n(i)+1:nm(i)],...
						[1:c(i) n(i)+1:nm(i)]);
	    for j=1:K
	      matrix_info.M{i,j}=matrix_info.M{i,j}([1:c(i) n(i)+1:nm(i)],...
						    [1:c(i) n(i)+1:nm(i)]); 
	    end;
        matrix_info.C{i} = matrix_info.C{i}(1:c(i),1:c(i));
        
	    n(i)=c(i);
	    matrix_info.n(i)=n(i);	    
	    nm(i)=c(i)+m(i);
	    matrix_info.nm(i)=nm(i);
        end;       
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do LQ-control if necessary or transform=1

for i= 1:N
        if n(i)~=0
        L{i}=zeros(m(i),n(i));
        if min(real(eig(matrix_info.A{i})))>=0|transform==1
	    disp('State-feedback performed')
            A=matrix_info.A{i};
	    B=matrix_info.B{i};
	    Q1=eye(n(i));
	    Q2=rho*eye(m(i));
	    [L{i},H,E]=lqr(A,B,Q1,Q2);
	    matrix_info.A{i}=matrix_info.A{i}-matrix_info.B{i}*L{i};
	    U{i}=[eye(n(i)) zeros(n(i),m(i));-L{i} eye(m(i))]; 
	    matrix_info.M0{i}=U{i}'*matrix_info.M0{i}*U{i};
            for j=1:K
	        matrix_info.M{i,j}=U{i}'*matrix_info.M{i,j}*U{i};
	    end;
        end; 
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check that the M matrices are linearly independent.
[checkM,matrix_info,Pbar,V]=check_M(matrix_info,solver,lowrank);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reformulate the problem with blockdiagonal A{i}.
% P = T^-T P_tilde T^-1

if lowrank
  T_diag = matrix_info.T_diag;
  T_diaginv = matrix_info.T_diaginv;
  D_diag = matrix_info.D_diag;
  
  for i = 1:N
    
    matrix_info.A{i} = D_diag{i};
    matrix_info.B{i} = T_diaginv{i}*matrix_info.B{i};
    matrix_info.C{i} = T_diaginv{i}*matrix_info.C{i}*T_diaginv{i}.';
    Konj_matrix = blkdiag(T_diag{i},eye(nm(i)-n(i)));
    
    matrix_info.M0{i} = ...
        Konj_matrix.'*matrix_info.M0{i}*Konj_matrix;
    
    for j = 1:K
      matrix_info.M{i,j} = ...
          Konj_matrix.'*matrix_info.M{i,j}*Konj_matrix;
    end
  end
end











