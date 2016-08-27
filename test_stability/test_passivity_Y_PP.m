%A B C D is given by admittance_to_state_matrix_PP.m
P = sdpvar(8,8);
F = [[P*A+A'*P+C'*C P*B+C'*D; B'*P+D'*C D'*D] >= 0];
F = [kyp(A,B,P) <= 0];
optimize(F);
Pfeasible = value(P);
checkset(F)
