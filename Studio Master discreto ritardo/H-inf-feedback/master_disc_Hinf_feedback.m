Mm = 0.64;
Bm = 12;

A = [0 1; 0 -Bm/Mm];
B = [0; 1/Mm];
C = [1 0; 0 1];
D = [0; 0];
B0 = [0; 1/Mm];
D0 = [0];

Ts = 0.001;
mu = 1;

csys = ss(A,B,C,D);
dsys = c2d(csys, Ts);

A = dsys.a;
B = dsys.b;
C = dsys.c;
D = dsys.d;

n = size(A,1);
p = size(B,1);
q = size(B,2);
r = size(D,1);
s = size(D,2);
t = size(C,1);
u = size(C,2);

P = sdpvar(n,n, 'symmetric');
H = sdpvar(n,n, 'symmetric');
L = sdpvar(q,n);
F = sdpvar(n,n); %
Q = sdpvar(n,n); %
R = sdpvar(q,t);
S = sdpvar(n,n); %
J = sdpvar(n,n); %
X = sdpvar(n,n);
Y = sdpvar(n,n); %

b11 = P;
b12 = J;
b13 = A*X+B*L;
b14 = A+B*R+C;
b15 = B+B*R*D;
b16 = zeros(n,n);

bb1 = [b11 b12 b13 b14 b15 b16];

b21 = b12';
b22 = H;
b23 = Q;
b24 = Y*A+F*C;
b25 = Y*B+F*D;
b26 = zeros(n,n);

bb2 = [b21 b22 b23 b24 b25 b26];

b31 = b13';
b32 = b23';
b33 = X+X'-P;
b34 = eye(n,n)+S'-J;
b35 = zeros(n,q);
b36 = X'*C'+L'*D';

bb3 = [b31 b32 b33 b34 b35 b36];

b41 = b14';
b42 = b24';
b43 = b34';
b44 = Y+Y'-H;
b45 = zeros(n,q);
b46 = C'+C'*R'*D';

bb4 = [b41 b42 b43 b44 b45 b46];

b51 = b15';
b52 = b25';
b53 = b35';
b54 = b45';
b55 = eye(q,q);
b56 = D'+D'*R'*D';

bb5 = [b51 b52 b53 b54 b55 b56];

b61 = b16';
b62 = b26';
b63 = b36';
b64 = b46';
b65 = b56';
b66 = mu*eye(n,n);

bb6 = [b61 b62 b63 b64 b65 b66];

LMI = [bb1; bb2; bb3; bb4; bb5; bb6];

res = [LMI >=0];

diagnostics = solvesdp(res);
disp(diagnostics.problem);
if diagnostics.problem == 0
 disp('Feasible')
 
 V = sqrt(-value(Y)*value(X)+value(S));
 U = V;
 m12 = -inv(V)*value(Y)*B;
 mat1 = [inv(V) m12; zeros(size(inv(V),1),size(inv(V),2)) eye(size(m12,1),size(m12,2))];
 mat2 = [value(Q)-value(Y)*A*value(X) value(F); value(L) value(R)];
 mat3 = [inv(value(U)) zeros(n,n); C*inv(value(X))*value(U) eye(n,n)];
 K = mat1*mat2*mat3;

elseif diagnostics.problem == 1
 disp('Infeasible')
else
 disp('Something else happened')
end
