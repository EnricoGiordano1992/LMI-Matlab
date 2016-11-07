Mm = 0.64;
Bm = 12;

A = [0 1; 0 -Bm/Mm];
B = [0; 1/Mm];
C = [1 0; 0 1];
D = [0; 0];
C0 = [1 10];
B0 = [0; 1/Mm];

n = size(A,1);
p = size(B,1);
q = size(B,2);
r = size(D,1);
s = size(D,2);

Ts = 0.01;

mu = 2;

csys = ss(A,B,C,D);
dsys = c2d(csys, Ts);

A = dsys.a;
B = dsys.b;
C = dsys.c;
D = dsys.d;

P = sdpvar(n,n, 'symmetric');
L = sdpvar(q,n);
X = sdpvar(n,n);

b11 = P;
b12 = A*X+B*L;
b13 = B;
b14 = zeros(n);
b21 = b12';
b22 = X+X'-P;
b23 = zeros(n, q);
b24 = X'*C'+L'*D';
b31 = b13';
b32 = b23';
b33 = eye(q, q);
b34 = D';
b41 = b14';
b42 = b24';
b43 = b34';
b44 = mu*eye(n,n);

lmi1 = [b11 b12 b13 b14];
lmi2 = [b21 b22 b23 b24];
lmi3 = [b31 b32 b33 b34];
lmi4 = [b41 b42 b43 b44];

LMI = [lmi1; lmi2; lmi3; lmi4];

F = [LMI >=0];

diagnostics = solvesdp(F);
disp(diagnostics.problem);
if diagnostics.problem == 0
 disp('Feasible')
 L_s = value(L);
 X_s = value(X);
 
 K = L_s * inv(X_s)

elseif diagnostics.problem == 1
 disp('Infeasible')
else
 disp('Something else happened')
end

