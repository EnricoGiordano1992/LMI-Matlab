function [U_diag,U,D,block_info] = create_blockdiag_sort(ua,ta)
% function [U_diag,D] = create_blockdiag(ua,ta)
% Blockdiagonalize a matrix A when eigenvalues and eigenvectors for A.' is known.
%
% [T,D] = eig(A') gives T and D
%
% Based  on the fact that A'*T = T*D
% <=>
% T'*A*inv(T)' = D
% U_diag = inv(T)'
%

n = size(ta,1);
D = zeros(n,n);%spalloc(n,n,2*n);
T_diag = zeros(n,n);%spalloc(n,n,n^2);
block_info = 0;

index = 1;
blk_index = 1;
real_index = n;
while index <= n
    eigenvalue = ta(index,index);
    if isreal(eigenvalue)
        D(real_index,real_index) = eigenvalue;
        T_diag(:,real_index) = ua(:,index);
        real_index = real_index-1;
        index = index+1;
    else
    a = real(eigenvalue);
    b = imag(eigenvalue);
    u = real(ua(:,index));
    v = imag(ua(:,index));
    D(blk_index:blk_index+1,blk_index:blk_index+1) = [a,-b;b,a];
    T_diag(:,blk_index) = u;
    T_diag(:,blk_index+1) = v;
    blk_index = blk_index + 2;
    index = index + 2;
    end
end
block_info = blk_index-1;

  U_diag = T_diag.'\eye(size(T_diag));
  U = T_diag.';