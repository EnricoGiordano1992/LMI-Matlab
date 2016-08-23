%A = [-1 2 0;-3 -4 1;0 0 -2];
%B = [1; 1; 1];
%C = [1 1 1];
%D = 0;
X = sdpvar(1,1);
F = [[[1 0; A B; C D]'*[0 X 0; X 0 0; 0 0 X]*[1 0; A B; C D]] >= 0];
F = [kyp(A,B,X) >= 0];
F = [F, trace(X) == 1];
optimize(F);
Pfeasible = value(X);
checkset(F)
