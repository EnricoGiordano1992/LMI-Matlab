%A B C D is given by admittance_to_state_matrix.m
P = sdpvar(1,1);
F = [[A'*P+P*A P*B-C'; B'*P-C -D'-D] > 0];
F = [kyp(A,B,P) < 0];
F = [F, trace(P) == 1];
optimize(F);
Pfeasible = value(P);
checkset(F)
