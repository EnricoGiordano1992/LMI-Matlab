Mm = 0.64;
Bm = 12;

s = tf('s');

A = [0 1; 0 -Bm/Mm];
B = [0; 1/Mm];
C = [1 0; 0 1];
D = [0; 0];
C0 = [1 10];
B0 = [0; 1/Mm];

n = size(A,1)
p = size(B,1)
q = size(B,2)

Q = sdpvar(n,n, 'symmetric');
M = sdpvar(q,p);
K = sdpvar(q,p);
E = zeros(q,q);
D0 = E;
Acl = A;

BlockUpLeft = Q*Acl'+M'*B'+Acl*Q+B*M;
BlockUpRight = B-Q*C';
BlockDownLeft = B'-C*Q;
BlockDownRight = zeros(size(BlockUpRight, 2), size(BlockDownLeft, 1));

LMI2 = [BlockUpLeft BlockUpRight; BlockDownLeft BlockDownRight];

F = [Q >= 0, LMI2 <=0];

diagnostics = solvesdp(F);
disp(diagnostics.problem);
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

