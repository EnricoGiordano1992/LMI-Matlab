n = size(A,1)
p = size(B,1)
q = size(B,2)
Q = sdpvar(n,n, 'symmetric');
M = sdpvar(q,p);

size(Q)
size(M)
size(A)
size(B)
size(C)

BlockUpLeft = Q*A'+M'*B'+A*Q+B*M;
BlockUpRight = B-Q*C';
BlockDownLeft = B'-C*Q;
BlockDownRight = zeros(size(BlockUpRight, 2), size(BlockDownLeft, 1));

F = [Q >= 0, [BlockUpLeft BlockUpRight; BlockDownLeft BlockDownRight] <=0]

diagnostics = solvesdp(F)
disp(diagnostics.problem)
if diagnostics.problem == 0
 disp('Feasible')
 Q_s = value(Q);
 M_s = value(M);
 
 K = M_s * Q_s'
 
elseif diagnostics.problem == 1
 disp('Infeasible')
else
 disp('Something else happened')
end
