function S=arraysum1(matrix_info,i,x)
  
K=matrix_info.K;
nm=matrix_info.nm(i);
M=matrix_info.M;

S=zeros(nm,nm);
for j=1:K
  S=S+x(j)*M{i,j};
end;


