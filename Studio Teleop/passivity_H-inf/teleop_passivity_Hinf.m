Ms = 0.64;
Bs = 12;
Mm = 0.64;
Bm = 12;

Ke = 0%0.01;
Me = 0%1;
De = 0%0.01; %damping
Xe = 9;

A_master = [0 1; 0 -Bm/Mm];
B_master = [0; 1/Mm];
C_master = [1 0; 0 1];
D_master = [0; 0];

A_slave = [0 1; 0 -Bs/Ms];
B_slave = [0; 1/Ms];
C_slave = [1 0; 0 1];
D_slave = [0; 0];

s = tf('s');

A = [0 1 0 0; 0 (-Bm/Mm) 0 0; 0 0 0 1; 0 0 0 -(Bs/Ms)];
B = [0 0; 1/Mm 0; 0 0; 0 1/Ms];
B0 = [0 0; 1/Mm 0; 0 0; 0 1/Ms];
%C0 = [1 0 -1 0; 0 1 0 -1];
%C0 = 4*[1 0 -1 0];
C0 = 4*[0 1 0 -1];
C = [0 1 0 0; 0 0 0 1];

n = size(A,1)
q = size(B,2)
r = size(C0,1)

Q = sdpvar(n,n, 'symmetric');
M = sdpvar(q,n);
K = sdpvar(q,n);
E = zeros(r,q);
D0 = E;
%Acl = A + B*K;
Acl = A;

gamma = 1;%0.000041;

b11 = (A*Q + B*M)'+(A*Q + B*M);
b12 = B0;
b13 = (C0*Q + E*M)';
b21 = B0';
b22 = -gamma^2 * eye(q);
b23 = D0';
b31 = (C0*Q + E*M);
b32 = D0;
b33 = -eye(r);

LMI1 = [b11 b12 b13; b21 b22 b23; b31 b32 b33];


BlockUpLeft = Q*Acl'+M'*B'+Acl*Q+B*M;
BlockUpRight = B-Q*C';
BlockDownLeft = B'-C*Q;
BlockDownRight = zeros(size(BlockUpRight, 2), size(BlockDownLeft, 1));

LMI2 = [BlockUpLeft BlockUpRight; BlockDownLeft BlockDownRight];

F = [Q >= 0, LMI1 <=0, LMI2 <= 0];

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

