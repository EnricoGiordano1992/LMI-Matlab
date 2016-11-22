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

gamma =2.2;%0.000041;

b11 = (A*Q + B*M)'+(A*Q + B*M);
b12 = B0;
b13 = (C0*Q + E*M)';
b13 = (C0*Q)';
b21 = B0';
b22 = -gamma^2 * eye(q);
b23 = D0';
b31 = (C0*Q + E*M);
b31 = (C0*Q);
b32 = D0;
b33 = -eye(q);

LMI1 = [b11 b12 b13; b21 b22 b23; b31 b32 b33];

F = [Q >= 0, LMI1 <=0];

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

