function [matrix_info,T] = transform_constraints(matrix_info_old)
% Transform the constraints with [T,0;0,I]  
  
  matrix_info = matrix_info_old;
  
  K = matrix_info.K;
  N = matrix_info.N;
  n = matrix_info.n;
  nm = matrix_info.nm;
  T = cell(N,1);
  for i = 1:N
    
    [T{i,1},D] = blockdiagonalize(matrix_info.A{i});
    
    matrix_info.A{i} = T{i,1}\matrix_info.A{i}*T{i,1};

    matrix_info.B{i} = T{i,1}\matrix_info.B{i};
    
    Konj_matrix = [T{i,1},zeros(n(i),nm(i)-n(i)); ...
                   zeros(nm(i)-n(i),n(i)),eye(nm(i)-n(i))];
    
    matrix_info.M0{i} = ...
        Konj_matrix.'*matrix_info.M0{i}*Konj_matrix;
                   
    for j = 1:K
      matrix_info.M{i,j} = ...
        Konj_matrix.'*matrix_info.M{i,j}*Konj_matrix;
    end
  end
  
                   