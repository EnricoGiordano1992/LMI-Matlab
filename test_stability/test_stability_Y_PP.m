%A is given by admittance_to_state_matrix_PP.m
P = sdpvar(6,6);

F = [P >= 0, A'*P+P*A <= 0];

F = [F, trace(P) == 1];

F

optimize(F);
Pfeasible = value(P);

checkset(F)

F = [P >= 0, A'*P+P*A <= 0, trace(P)==1];
optimize(F,P(1,1));

