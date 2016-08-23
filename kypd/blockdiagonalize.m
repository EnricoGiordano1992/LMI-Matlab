function [T,D] = blockdiagonalize(matrix)
% Gives the transformationmatrices T and D that are defined by:
% 
% T^-1*A*T = D 
%  
% where D is a block-diagonal matrix with real and 2x2
% matrices in the diagonal. The singular values are real
% eigenvalues to A and the matrices are connected with complex
% eigenvalues as eigvalues = a +- b*i <=> [a,b;-b,a] in the diagonal.
% Corresponding eigenvector are u +- v*i
  
  n = size(matrix,1);
  
  [T,D] = eig(matrix,eye(n),'qz');
  
  i=1;
  while i <= n
    if ~isreal(D(i,i))
      a = real(D(i,i));
      b = imag(D(i,i));
      u = real(T(:,i));
      v = imag(T(:,i));
      D(i:i+1,i:i+1) = [a,b;-b,a];
      T(:,i) = u;
      T(:,i+1) = v;
      i = i+2;
    else     
      i = i+1;
    end
  end
      