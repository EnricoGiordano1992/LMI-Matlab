function schurpar = derankcross_sort(Basis_matrices,len,block_end,n,nm)
rcum = ones(len+1,1);
alphas = [];
V = [];

%Make 'derankcross' for double-crosses
 for i = 1:(nm-n)*block_end %number of dubble-crosses
     Fi = Basis_matrices{i};
     block_number = ceil(i/(2*(nm-n)));
     non_zero = [2*block_number-1 2*block_number];
     Cross1 = zeros(size(Fi));
     Cross2 = Cross1;
     Cross1(non_zero(1),:) = Fi(non_zero(1),:);
     Cross1(:,non_zero(1)) = Fi(:,non_zero(1));
     Cross2(non_zero(2),:) = Fi(non_zero(2),:);
     Cross2(:,non_zero(2)) = Fi(:,non_zero(2));
     Cross1(non_zero(1),non_zero(2)) = Cross1(non_zero(1),non_zero(2))*.5;
     Cross1(non_zero(2),non_zero(1)) = Cross1(non_zero(2),non_zero(1))*.5;
     Cross2(non_zero(1),non_zero(2)) = Cross2(non_zero(1),non_zero(2))*.5;
     Cross2(non_zero(2),non_zero(1)) = Cross2(non_zero(2),non_zero(1))*.5;        
     [alpha_temp V_temp] = derankcross_single(sparse(Cross1),non_zero(1));
     alphas = [alphas, alpha_temp];
     V = [V, V_temp];
     [alpha_temp V_temp] = derankcross_single(sparse(Cross2),non_zero(2));
     alphas = [alphas, alpha_temp];
     V = [V, V_temp];
     rcum(i+1) = rcum(i)+4;
 end
 
 %Single-crosses
 for i = 1:(nm-n)*(n-block_end)%number of single-crosses
     position = (nm-n)*block_end+i;
     Fi = Basis_matrices{position};
     Cross = zeros(size(Fi));
     non_zero = block_end+ceil(i/(nm-n));
     Cross(:,non_zero) = Fi(:,non_zero);
     Cross(non_zero,:) = Fi(non_zero,:);
     [alpha_temp V_temp] = derankcross_single(sparse(Cross),non_zero);
     rcum(position+1) = rcum(position) + 2;
     alphas = [alphas, alpha_temp];
     V = [V, V_temp];
 end
 
 %The 'extra block'
 i = (nm-n)*n+1;
 for l=1:nm-n
     for k=l:nm-n
         if k==l
             rcum(i+1) = rcum(i) + 1;
             alphas = [alphas, 1];
             v = zeros(nm,1);
             v(n+k) = 1;
             V = [V,v];
             i=i+1;
         else
             Fi = Basis_matrices{i};
             [alpha_temp V_temp] = derankcross_single(sparse(Fi),n+k);
             rcum(i+1) = rcum(i) + 2;
             alphas = [alphas, alpha_temp];
             V = [V, V_temp];
             i=i+1;
         end;
     end;
 end;

schurpar = {{rcum alphas V}};





function [alpha V] = derankcross_single(Cross, diag)
%function [alpha V] = derankcross_single(Cross, diag)
% 
% Rewrite a cross-matrix with rank 1 as alpha*V*V.'
%
% Cross is a matrix with non-zero element in position located on row/column
% 'diag'.
%

v = zeros(size(Cross,1),1);
v(:,1) = Cross(:,diag);
v(diag) = .5*v(diag);
normv = norm(v);
v1 = 1/normv*v;
v2 = v1;
v1(diag) = v1(diag)+1;
v2(diag) = v2(diag)-1;
alpha = [.5*normv -.5*normv];
V = [v1,v2];

