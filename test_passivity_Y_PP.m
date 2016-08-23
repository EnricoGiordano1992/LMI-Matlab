%A B C D is given by admittance_to_state_matrix_PP.m
P = sdpvar(6,6);
F = [[P*A+A'*P+C'*C P*B+C'*D; B'*P+D'*C D'*D] >= 0];
F = [kyp(A,B,P) <= 0];
F = [F, trace(P) == 1];
optimize(F);
Pfeasible = value(P);
checkset(F)
