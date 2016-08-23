function [F,n_basis,ntot_basis,F0]=basis_matrices(matrix_info,solver,lowrank);

% Calculates all basis matrices for Z.
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2|isempty(solver)
     solver='schur';
end;

N=matrix_info.N;
n=matrix_info.n;
nm=matrix_info.nm;
A=matrix_info.A;
B=matrix_info.B;
F0 = cell(N,1);
t=1;
for i=1:N
    F{t}=spalloc(nm(i),nm(i),nm(i)*nm(i));
    temp=t;
    if n(i)~=0
      A1=A{i};
      [ma,na] = size(A1);
      [mb,nb] = size(A1');
      if solver=='schur'
          [ua,ta] = schur(A1,'complex');
          j = ma:-1:1;
          ub = ua(:,j);
          tb = ta(j,j)';
      end;
      
      if solver=='diago'
          [ua,ta] = eig(A1);
          uainv = inv(ua);
      end;
    end;

    % Solve lyap(F0{i}) = Ci

    C = matrix_info.C{i};
    F0{i} = spalloc(nm(i),nm(i),nm(i)*nm(i));

    if solver=='schur'
        F0{i}(1:n(i),1:n(i)) = lyapunov(ua,ta,ub,tb,-C);
    end;

    if solver=='diago'
        F0{i}(1:n(i),1:n(i)) = diaglyapunov(ua,uainv,ta,-C);
    end;


    for k=1:n(i) % number of rows i B (A)
        for l=1:nm(i)-n(i) % number of columns i B
      
        F{t}=spalloc(nm(i),nm(i),nm(i)*nm(i));
        F{t}(n(i)+l,k)=1;
        F{t}(k,n(i)+l)=1;
        
        if solver=='schur'
          F{t}(1:n(i),1:n(i))=lyapunov(ua,ta,ub,tb,B{i}*F{t}(1:n(i),...
                                    n(i)+1:nm(i))'+...
                                    F{t}(1:n(i),n(i)+ 1:nm(i))* ...
				       B{i}');
        end;
        
        if solver=='diago'
            F{t}(1:n(i),1:n(i))=diaglyapunov(ua,uainv,ta,B{i}*F{t}(1:n(i),...
                                    n(i)+1:nm(i))'+...
                                    F{t}(1:n(i),n(i)+ 1:nm(i))* ...
                                             B{i}');
        end;
        t=t+1;
      end;
    end;
    for l=1:nm(i)-n(i)
        for k=l:nm(i)-n(i)
            F{t}=spalloc(nm(i),nm(i),nm(i)*nm(i));
            F{t}(n(i)+l,n(i)+k)=1;
            if k~=l
                F{t}(n(i)+k,n(i)+l)=1;
            end;
            t=t+1;
        end;
    end;
    n_basis(i)=t-temp;
end;

    
ntot_basis=sum(n_basis);
